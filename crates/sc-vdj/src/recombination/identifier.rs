use super::Recombination;
use crate::index::{Chain, SegmentKind, VdjIndex};
use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReceptorRole {
    Heavy,
    Light,
}
impl ReceptorRole {
    fn prefix(self) -> &'static str {
        if self == Self::Heavy {
            "HC"
        } else {
            "LC"
        }
    }
    fn from_prefix(s: &str) -> Option<Self> {
        match s {
            "HC" => Some(Self::Heavy),
            "LC" => Some(Self::Light),
            _ => None,
        }
    }
}

/// Compact nucleotide-independent structural recombination identity.
///
/// Version 3 keeps the structural identity semantics while making the segment
/// part self-describing. V/D/J are global `SegmentId`s from the VDJ index's
/// shared segment pool. One width nibble says how many hex digits each ID uses,
/// so v3 can always be decoded numerically without loading the index.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RecombinationId {
    role: ReceptorRole,
    payload: String,
}
impl RecombinationId {
    const VERSION: u8 = 3;
    const LEGACY_V2: u8 = 2;
    pub(crate) fn placeholder(chain: Chain) -> Self {
        Self {
            role: if chain.has_d() {
                ReceptorRole::Heavy
            } else {
                ReceptorRole::Light
            },
            payload: String::new(),
        }
    }
    pub fn from_recombination(r: &Recombination, _index: &VdjIndex) -> Result<Self, String> {
        let role = if r.chain.has_d() {
            ReceptorRole::Heavy
        } else {
            ReceptorRole::Light
        };
        let d = if role == ReceptorRole::Heavy {
            Some(r.d.ok_or("heavy-chain recombination lacks D")?)
        } else {
            None
        };
        let width = global_id_width(r.v, d, r.j);
        let mut p = format!("{:X}{:X}{:X}", Self::VERSION, r.chain.code(), width);
        push_global(&mut p, r.v, width);
        if let Some(d) = d {
            push_global(&mut p, d, width);
        }
        push_global(&mut p, r.j, width);
        let j = &r.junction;
        push_measurement(&mut p, j.v_del_3);
        push_measurement(&mut p, j.p_v3_len());
        push_measurement(&mut p, j.n1_len());
        if role == ReceptorRole::Heavy {
            push_measurement(&mut p, j.p_d5_len());
            push_measurement(&mut p, j.d_del_5.ok_or("missing d_del_5")?);
            push_measurement(&mut p, j.d_retained_len());
            push_measurement(&mut p, j.d_del_3.ok_or("missing d_del_3")?);
            push_measurement(&mut p, j.p_d3_len());
            push_measurement(&mut p, j.n2_len());
        }
        push_measurement(&mut p, j.p_j5_len());
        push_measurement(&mut p, j.j_del_5);
        p.push(if j.pn_alternative { '1' } else { '0' });
        Ok(Self { role, payload: p })
    }
    /// Decode the structural fields that are intrinsic to the compact ID.
    ///
    /// Version 1 and version 3 IDs decode uniquely without a reference index.
    /// Legacy version 2 stored variable-width chain/kind-local ordinals; without
    /// the index that supplied those widths, more than one V/D/J split can be
    /// valid. In that case all valid legacy candidates are returned rather than
    /// guessed.
    pub fn decode_numeric_candidates(&self) -> Result<Vec<DecodedNumericRecombinationId>, String> {
        let mut c = Cursor::new(&self.payload);
        let version = c.nibble()?;
        match version {
            1 => self.decode_v1_numeric(c).map(|x| vec![x]),
            Self::LEGACY_V2 => self.decode_v2_numeric_candidates(c),
            Self::VERSION => self.decode_v3_numeric(c).map(|x| vec![x]),
            _ => Err(format!(
                "unsupported compact recombination id version {version}"
            )),
        }
    }

    fn decode_v1_numeric(
        &self,
        mut c: Cursor<'_>,
    ) -> Result<DecodedNumericRecombinationId, String> {
        let v_id = c.hex(3)?;
        let d_id = if self.role == ReceptorRole::Heavy {
            Some(c.hex(3)?)
        } else {
            None
        };
        let j_id = c.hex(3)?;
        finish_numeric_decode(&mut c, None, v_id, d_id, j_id)
    }

    fn decode_v3_numeric(
        &self,
        mut c: Cursor<'_>,
    ) -> Result<DecodedNumericRecombinationId, String> {
        let chain = Chain::from_code(c.nibble()?).map_err(|e| e.to_string())?;
        if chain.has_d() != (self.role == ReceptorRole::Heavy) {
            return Err(format!("{} prefix incompatible with {chain}", self.role.prefix()));
        }
        let width = c.nibble()? as usize;
        if !(1..=4).contains(&width) {
            return Err(format!("invalid v3 segment-id width {width}"));
        }
        let v_id = c.hex(width)?;
        let d_id = if chain.has_d() { Some(c.hex(width)?) } else { None };
        let j_id = c.hex(width)?;
        finish_numeric_decode(&mut c, Some(chain), v_id, d_id, j_id)
    }

    fn decode_v2_numeric_candidates(
        &self,
        c: Cursor<'_>,
    ) -> Result<Vec<DecodedNumericRecombinationId>, String> {
        let mut header = c;
        let chain = Chain::from_code(header.nibble()?).map_err(|e| e.to_string())?;
        if chain.has_d() != (self.role == ReceptorRole::Heavy) {
            return Err(format!(
                "{} prefix incompatible with {chain}",
                self.role.prefix()
            ));
        }

        let start = header.p;
        let d_widths: &[usize] = if chain.has_d() { &[1, 2, 3, 4] } else { &[0] };
        let mut out = Vec::new();
        for v_width in 1..=4 {
            for &d_width in d_widths {
                for j_width in 1..=4 {
                    let mut candidate = Cursor {
                        b: header.b,
                        p: start,
                    };
                    let parsed = (|| {
                        let v_id = candidate.hex(v_width)?;
                        let d_id = if chain.has_d() {
                            Some(candidate.hex(d_width)?)
                        } else {
                            None
                        };
                        let j_id = candidate.hex(j_width)?;
                        finish_numeric_decode(&mut candidate, Some(chain), v_id, d_id, j_id)
                    })();
                    if let Ok(decoded) = parsed {
                        if !out.contains(&decoded) {
                            out.push(decoded);
                        }
                    }
                }
            }
        }
        if out.is_empty() {
            return Err("no valid index-free V/D/J split for compact recombination id".into());
        }
        Ok(out)
    }

    pub fn decode(&self, index: &VdjIndex) -> Result<DecodedRecombinationId, String> {
        let mut c = Cursor::new(&self.payload);
        let version = c.nibble()?;
        match version {
            1 => self.decode_v1(index, c),
            Self::LEGACY_V2 => self.decode_v2(index, c),
            Self::VERSION => self.decode_v3(index, c),
            _ => Err(format!(
                "unsupported compact recombination id version {version}"
            )),
        }
    }
    fn decode_v3(
        &self,
        index: &VdjIndex,
        mut c: Cursor<'_>,
    ) -> Result<DecodedRecombinationId, String> {
        let chain = Chain::from_code(c.nibble()?).map_err(|e| e.to_string())?;
        if chain.has_d() != (self.role == ReceptorRole::Heavy) {
            return Err(format!("{} prefix incompatible with {chain}", self.role.prefix()));
        }
        let width = c.nibble()? as usize;
        if !(1..=4).contains(&width) {
            return Err(format!("invalid v3 segment-id width {width}"));
        }
        let v_id = c.hex(width)?;
        let d_id = if chain.has_d() { Some(c.hex(width)?) } else { None };
        let j_id = c.hex(width)?;
        let v = checked_global(index, v_id, chain, SegmentKind::V)?;
        let d = d_id
            .map(|id| checked_global(index, id, chain, SegmentKind::D))
            .transpose()?;
        let j = checked_global(index, j_id, chain, SegmentKind::J)?;
        finish_decode(self.role, index, &mut c, chain, v, d, j)
    }

    fn decode_v2(
        &self,
        index: &VdjIndex,
        mut c: Cursor<'_>,
    ) -> Result<DecodedRecombinationId, String> {
        let chain = Chain::from_code(c.nibble()?).map_err(|e| e.to_string())?;
        if chain.has_d() != (self.role == ReceptorRole::Heavy) {
            return Err(format!(
                "{} prefix incompatible with {chain}",
                self.role.prefix()
            ));
        }
        let v = read_local(&mut c, index, chain, SegmentKind::V)?;
        let d = if chain.has_d() {
            Some(read_local(&mut c, index, chain, SegmentKind::D)?)
        } else {
            None
        };
        let j = read_local(&mut c, index, chain, SegmentKind::J)?;
        finish_decode(self.role, index, &mut c, chain, v, d, j)
    }
    fn decode_v1(
        &self,
        index: &VdjIndex,
        mut c: Cursor<'_>,
    ) -> Result<DecodedRecombinationId, String> {
        let v_index = c.hex(3)?;
        let d_index = if self.role == ReceptorRole::Heavy {
            Some(c.hex(3)?)
        } else {
            None
        };
        let j_index = c.hex(3)?;
        let v = index
            .segments
            .get(v_index)
            .ok_or_else(|| format!("V segment index {v_index} outside index"))?;
        let j = index
            .segments
            .get(j_index)
            .ok_or_else(|| format!("J segment index {j_index} outside index"))?;
        if v.kind != SegmentKind::V || j.kind != SegmentKind::J || v.chain != j.chain {
            return Err("v1 V/J indices do not resolve to one valid chain".into());
        }
        let chain = v.chain;
        if chain.has_d() != (self.role == ReceptorRole::Heavy) {
            return Err(format!(
                "{} prefix incompatible with {chain}",
                self.role.prefix()
            ));
        }
        let d = if let Some(i) = d_index {
            let x = index
                .segments
                .get(i)
                .ok_or_else(|| format!("D segment index {i} outside index"))?;
            if x.kind != SegmentKind::D || x.chain != chain {
                return Err("v1 D index is not a matching D segment".into());
            }
            Some(x)
        } else {
            None
        };
        finish_decode(self.role, index, &mut c, chain, v, d, j)
    }
    /// Re-encode a legacy v1/v2 identifier as the self-describing v3 format.
    /// The index is only required to translate legacy segment ordinals/IDs.
    pub fn to_v3(&self, index: &VdjIndex) -> Result<Self, String> {
        let mut c = Cursor::new(&self.payload);
        let version = c.nibble()?;
        if version == Self::VERSION {
            self.decode(index)?;
            return Ok(self.clone());
        }

        let numeric = match version {
            1 => self.decode_v1_global_numeric(index, c)?,
            Self::LEGACY_V2 => self.decode_v2_global_numeric(index, c)?,
            _ => return Err(format!("unsupported compact recombination id version {version}")),
        };
        Self::from_global_numeric(self.role, &numeric)
    }

    fn decode_v1_global_numeric(
        &self,
        index: &VdjIndex,
        mut c: Cursor<'_>,
    ) -> Result<DecodedNumericRecombinationId, String> {
        let v_id = c.hex(3)?;
        let d_id = if self.role == ReceptorRole::Heavy { Some(c.hex(3)?) } else { None };
        let j_id = c.hex(3)?;
        let v = checked_global(index, v_id, index.segment(v_id as u16).ok_or("V segment outside index")?.chain, SegmentKind::V)?;
        let chain = v.chain;
        checked_global(index, j_id, chain, SegmentKind::J)?;
        if let Some(id) = d_id { checked_global(index, id, chain, SegmentKind::D)?; }
        finish_numeric_decode(&mut c, Some(chain), v_id, d_id, j_id)
    }

    fn decode_v2_global_numeric(
        &self,
        index: &VdjIndex,
        mut c: Cursor<'_>,
    ) -> Result<DecodedNumericRecombinationId, String> {
        let chain = Chain::from_code(c.nibble()?).map_err(|e| e.to_string())?;
        if chain.has_d() != (self.role == ReceptorRole::Heavy) {
            return Err(format!("{} prefix incompatible with {chain}", self.role.prefix()));
        }
        let v = read_local(&mut c, index, chain, SegmentKind::V)?.id as usize;
        let d = if chain.has_d() {
            Some(read_local(&mut c, index, chain, SegmentKind::D)?.id as usize)
        } else {
            None
        };
        let j = read_local(&mut c, index, chain, SegmentKind::J)?.id as usize;
        finish_numeric_decode(&mut c, Some(chain), v, d, j)
    }

    fn from_global_numeric(
        role: ReceptorRole,
        d: &DecodedNumericRecombinationId,
    ) -> Result<Self, String> {
        let chain = d.chain.ok_or("legacy identifier does not contain a chain")?;
        let v = u16::try_from(d.v_id).map_err(|_| "V segment ID exceeds u16")?;
        let dd = d.d_id.map(u16::try_from).transpose().map_err(|_| "D segment ID exceeds u16")?;
        let j = u16::try_from(d.j_id).map_err(|_| "J segment ID exceeds u16")?;
        let width = global_id_width(v, dd, j);
        let mut p = format!("{:X}{:X}{:X}", Self::VERSION, chain.code(), width);
        push_global(&mut p, v, width);
        if let Some(dd) = dd { push_global(&mut p, dd, width); }
        push_global(&mut p, j, width);
        push_measurement(&mut p, d.v_del_3);
        push_measurement(&mut p, d.p_v3_len);
        push_measurement(&mut p, d.n1_len);
        if chain.has_d() {
            push_measurement(&mut p, d.p_d5_len.ok_or("missing p_d5_len")?);
            push_measurement(&mut p, d.d_del_5.ok_or("missing d_del_5")?);
            push_measurement(&mut p, d.d_retained_len.ok_or("missing d_retained_len")?);
            push_measurement(&mut p, d.d_del_3.ok_or("missing d_del_3")?);
            push_measurement(&mut p, d.p_d3_len.ok_or("missing p_d3_len")?);
            push_measurement(&mut p, d.n2_len.ok_or("missing n2_len")?);
        }
        push_measurement(&mut p, d.p_j5_len);
        push_measurement(&mut p, d.j_del_5);
        p.push(if d.pn_alternative { '1' } else { '0' });
        Ok(Self { role, payload: p })
    }

    pub fn payload(&self) -> &str {
        &self.payload
    }
}
impl fmt::Display for RecombinationId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.role.prefix(), self.payload)
    }
}
impl FromStr for RecombinationId {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        let (p, x) = s
            .split_once(':')
            .ok_or("recombination id requires HC:/LC: prefix")?;
        let role = ReceptorRole::from_prefix(p).ok_or_else(|| format!("unknown prefix {p}"))?;
        if x.is_empty() || !x.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err("payload must be hexadecimal".into());
        }
        Ok(Self {
            role,
            payload: x.to_ascii_uppercase(),
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedNumericRecombinationId {
    /// Chain is encoded directly in version 2 IDs. Version 1 needs an index to
    /// recover the chain, so this is `None` for index-free v1 decoding.
    pub chain: Option<Chain>,
    /// Version 2: chain-local V ordinal. Version 1: global segment index.
    pub v_id: usize,
    /// Version 2: chain-local D ordinal. Version 1: global segment index.
    pub d_id: Option<usize>,
    /// Version 2: chain-local J ordinal. Version 1: global segment index.
    pub j_id: usize,
    pub v_del_3: u16,
    pub p_v3_len: u16,
    pub n1_len: u16,
    pub p_d5_len: Option<u16>,
    pub d_del_5: Option<u16>,
    pub d_retained_len: Option<u16>,
    pub d_del_3: Option<u16>,
    pub p_d3_len: Option<u16>,
    pub n2_len: Option<u16>,
    pub p_j5_len: u16,
    pub j_del_5: u16,
    pub pn_alternative: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedRecombinationId {
    pub chain: Chain,
    pub v: String,
    pub d: Option<String>,
    pub j: String,
    pub v_del_3: u16,
    pub p_v3_len: u16,
    pub n1_len: u16,
    pub p_d5_len: Option<u16>,
    pub d_del_5: Option<u16>,
    pub d_retained_len: Option<u16>,
    pub d_del_3: Option<u16>,
    pub p_d3_len: Option<u16>,
    pub n2_len: Option<u16>,
    pub p_j5_len: u16,
    pub j_del_5: u16,
    pub pn_alternative: bool,
}

fn finish_numeric_decode(
    c: &mut Cursor<'_>,
    chain: Option<Chain>,
    v_id: usize,
    d_id: Option<usize>,
    j_id: usize,
) -> Result<DecodedNumericRecombinationId, String> {
    let has_d = d_id.is_some();
    let v_del_3 = c.measurement()?;
    let p_v3_len = c.measurement()?;
    let n1_len = c.measurement()?;
    let (p_d5_len, d_del_5, d_retained_len, d_del_3, p_d3_len, n2_len) = if has_d {
        (
            Some(c.measurement()?),
            Some(c.measurement()?),
            Some(c.measurement()?),
            Some(c.measurement()?),
            Some(c.measurement()?),
            Some(c.measurement()?),
        )
    } else {
        (None, None, None, None, None, None)
    };
    let p_j5_len = c.measurement()?;
    let j_del_5 = c.measurement()?;
    let pn_alternative = match c.nibble()? {
        0 => false,
        1 => true,
        x => return Err(format!("invalid P/N flag {x}")),
    };
    if !c.finished() {
        return Err("trailing hex data".into());
    }
    Ok(DecodedNumericRecombinationId {
        chain,
        v_id,
        d_id,
        j_id,
        v_del_3,
        p_v3_len,
        n1_len,
        p_d5_len,
        d_del_5,
        d_retained_len,
        d_del_3,
        p_d3_len,
        n2_len,
        p_j5_len,
        j_del_5,
        pn_alternative,
    })
}

fn finish_decode(
    _role: ReceptorRole,
    _index: &VdjIndex,
    c: &mut Cursor<'_>,
    chain: Chain,
    v: &crate::index::VdjSegment,
    d: Option<&crate::index::VdjSegment>,
    j: &crate::index::VdjSegment,
) -> Result<DecodedRecombinationId, String> {
    let v_del_3 = c.measurement()?;
    let p_v3_len = c.measurement()?;
    let n1_len = c.measurement()?;
    let (p_d5_len, d_del_5, d_retained_len, d_del_3, p_d3_len, n2_len) = if chain.has_d() {
        (
            Some(c.measurement()?),
            Some(c.measurement()?),
            Some(c.measurement()?),
            Some(c.measurement()?),
            Some(c.measurement()?),
            Some(c.measurement()?),
        )
    } else {
        (None, None, None, None, None, None)
    };
    let p_j5_len = c.measurement()?;
    let j_del_5 = c.measurement()?;
    let pn_alternative = match c.nibble()? {
        0 => false,
        1 => true,
        x => return Err(format!("invalid P/N flag {x}")),
    };
    if !c.finished() {
        return Err("trailing hex data".into());
    }
    Ok(DecodedRecombinationId {
        chain,
        v: v.name.clone(),
        d: d.map(|x| x.name.clone()),
        j: j.name.clone(),
        v_del_3,
        p_v3_len,
        n1_len,
        p_d5_len,
        d_del_5,
        d_retained_len,
        d_del_3,
        p_d3_len,
        n2_len,
        p_j5_len,
        j_del_5,
        pn_alternative,
    })
}

fn global_id_width(v: u16, d: Option<u16>, j: u16) -> usize {
    let max_id = d.map_or(v.max(j), |d| v.max(d).max(j));
    if max_id <= 0xF { 1 } else if max_id <= 0xFF { 2 } else if max_id <= 0xFFF { 3 } else { 4 }
}

fn push_global(out: &mut String, id: u16, width: usize) {
    out.push_str(&format!("{id:0width$X}"));
}

fn checked_global<'a>(
    index: &'a VdjIndex,
    id: usize,
    chain: Chain,
    kind: SegmentKind,
) -> Result<&'a crate::index::VdjSegment, String> {
    let id = u16::try_from(id).map_err(|_| format!("segment ID {id} exceeds u16"))?;
    let segment = index.segment(id).ok_or_else(|| format!("segment ID {id} outside index"))?;
    if segment.chain != chain || segment.kind != kind {
        return Err(format!("segment ID {id} is {} {:?}, expected {chain} {kind:?}", segment.chain, segment.kind));
    }
    Ok(segment)
}

fn hex_width(count: usize) -> usize {
    let mut n = count.saturating_sub(1);
    let mut w = 1;
    while n >= 16 {
        n >>= 4;
        w += 1
    }
    w
}
fn read_local<'a>(
    c: &mut Cursor<'_>,
    index: &'a VdjIndex,
    chain: Chain,
    kind: SegmentKind,
) -> Result<&'a crate::index::VdjSegment, String> {
    let w = hex_width(index.local_count(chain, kind));
    let ord = c.hex(w)?;
    index
        .segment_by_local_ordinal(chain, kind, ord)
        .ok_or_else(|| format!("{chain} {kind:?} local ordinal {ord} outside index"))
}
fn push_measurement(out: &mut String, v: u16) {
    if v <= 14 {
        out.push(char::from_digit(v as u32, 16).unwrap().to_ascii_uppercase())
    } else {
        out.push('F');
        out.push_str(&format!("{v:04X}"))
    }
}
struct Cursor<'a> {
    b: &'a [u8],
    p: usize,
}
impl<'a> Cursor<'a> {
    fn new(s: &'a str) -> Self {
        Self {
            b: s.as_bytes(),
            p: 0,
        }
    }
    fn nibble(&mut self) -> Result<u8, String> {
        let x = *self.b.get(self.p).ok_or("truncated recombination id")?;
        self.p += 1;
        match x {
            b'0'..=b'9' => Ok(x - b'0'),
            b'A'..=b'F' => Ok(x - b'A' + 10),
            b'a'..=b'f' => Ok(x - b'a' + 10),
            _ => Err("invalid hex".into()),
        }
    }
    fn hex(&mut self, n: usize) -> Result<usize, String> {
        let mut x = 0;
        for _ in 0..n {
            x = (x << 4) | self.nibble()? as usize;
        }
        Ok(x)
    }
    fn measurement(&mut self) -> Result<u16, String> {
        let x = self.nibble()?;
        if x < 15 {
            return Ok(x as u16);
        }
        let mut v = 0;
        for _ in 0..4 {
            v = (v << 4) | self.nibble()? as u16
        }
        Ok(v)
    }
    fn finished(&self) -> bool {
        self.p == self.b.len()
    }
}
