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
/// Version 2 keeps the original identity semantics, but V/D/J identifiers are
/// local ordinals inside the already-known chain/kind table.  The index tells
/// the decoder exactly how many hex digits each local ordinal requires.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RecombinationId {
    role: ReceptorRole,
    payload: String,
}
impl RecombinationId {
    const VERSION: u8 = 2;
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
    pub fn from_recombination(r: &Recombination, index: &VdjIndex) -> Result<Self, String> {
        let role = if r.chain.has_d() {
            ReceptorRole::Heavy
        } else {
            ReceptorRole::Light
        };
        let mut p = format!("{:X}{:X}", Self::VERSION, r.chain.code());
        push_local(&mut p, index, r.chain, SegmentKind::V, r.v)?;
        if role == ReceptorRole::Heavy {
            push_local(
                &mut p,
                index,
                r.chain,
                SegmentKind::D,
                r.d.ok_or("heavy-chain recombination lacks D")?,
            )?;
        }
        push_local(&mut p, index, r.chain, SegmentKind::J, r.j)?;
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
    pub fn decode(&self, index: &VdjIndex) -> Result<DecodedRecombinationId, String> {
        let mut c = Cursor::new(&self.payload);
        let version = c.nibble()?;
        match version {
            1 => self.decode_v1(index, c),
            Self::VERSION => self.decode_v2(index, c),
            _ => Err(format!(
                "unsupported compact recombination id version {version}"
            )),
        }
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

fn hex_width(count: usize) -> usize {
    let mut n = count.saturating_sub(1);
    let mut w = 1;
    while n >= 16 {
        n >>= 4;
        w += 1
    }
    w
}
fn push_local(
    out: &mut String,
    index: &VdjIndex,
    chain: Chain,
    kind: SegmentKind,
    id: u16,
) -> Result<(), String> {
    let ord = index
        .local_ordinal(id)
        .ok_or_else(|| format!("segment {id} not present in local table"))?;
    let w = hex_width(index.local_count(chain, kind));
    if ord >= 1usize << (4 * w) {
        return Err("local segment ordinal exceeds encoded width".into());
    }
    out.push_str(&format!("{ord:0w$X}"));
    Ok(())
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
