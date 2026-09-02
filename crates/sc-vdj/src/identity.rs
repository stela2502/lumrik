use crate::posterior::{CellVdjSummary, RearrangementCall};
use crate::types::Chain;
use std::fmt;
use std::str::FromStr;

/// Nucleotide-independent structural description of one V(D)J recombination.
///
/// The selected germline segment names and junction geometry are retained, but
/// the actual N/P/junction nucleotide sequence is deliberately not part of this
/// identity. `None` means that a measurement has not yet been resolved; fields
/// that do not apply to VJ chains are omitted from the canonical code.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RecombinationMeasurements {
    pub chain: Chain,
    pub v: String,
    pub v_del_3: Option<u16>,
    pub p_v3_len: Option<u16>,
    pub n1_len: Option<u16>,
    pub p_d5_len: Option<u16>,
    pub d: Option<String>,
    pub d_del_5: Option<u16>,
    pub d_retained_len: Option<u16>,
    pub d_del_3: Option<u16>,
    pub p_d3_len: Option<u16>,
    pub n2_len: Option<u16>,
    pub p_j5_len: Option<u16>,
    pub j_del_5: Option<u16>,
    pub j: String,
    pub pn_alternative: bool,
}

impl RecombinationMeasurements {
    /// Create the structural identity currently available from a posterior call.
    /// Measured junction geometry is used when one continuous coherent UMI contig
    /// supports the selected segments; unresolved fields stay explicit as `?`.
    pub fn from_call(call: &RearrangementCall) -> Option<Self> {
        let v = call.v.as_ref()?.id.clone();
        let j = call.j.as_ref()?.id.clone();
        let junction = call.junction.as_ref();
        Some(Self {
            chain: call.chain,
            v,
            v_del_3: junction.map(|x| x.v_del_3),
            p_v3_len: junction.map(|x| x.p_v3_len()),
            n1_len: junction.map(|x| x.n1_len()),
            p_d5_len: if call.chain.has_d() { junction.map(|x| x.p_d5_len()) } else { None },
            d: call.d.as_ref().map(|x| x.id.clone()),
            d_del_5: junction.and_then(|x| x.d_del_5),
            d_retained_len: junction.and_then(|x| x.d_retained_len),
            d_del_3: junction.and_then(|x| x.d_del_3),
            p_d3_len: if call.chain.has_d() { junction.map(|x| x.p_d3_len()) } else { None },
            n2_len: if call.chain.has_d() { junction.map(|x| x.n2_len()) } else { None },
            p_j5_len: junction.map(|x| x.p_j5_len()),
            j_del_5: junction.map(|x| x.j_del_5),
            j,
            pn_alternative: junction.is_some_and(|x| x.pn_alternative),
        })
    }


    pub fn is_complete(&self) -> bool {
        let common = self.v_del_3.is_some()
            && self.p_v3_len.is_some()
            && self.n1_len.is_some()
            && self.p_j5_len.is_some()
            && self.j_del_5.is_some();
        if !common {
            return false;
        }
        if self.chain.has_d() {
            self.p_d5_len.is_some()
                && self.d_del_5.is_some()
                && self.d_retained_len.is_some()
                && self.d_del_3.is_some()
                && self.p_d3_len.is_some()
                && self.n2_len.is_some()
        } else {
            true
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReceptorRole {
    Heavy,
    Light,
}

impl ReceptorRole {
    pub fn prefix(self) -> &'static str {
        match self {
            Self::Heavy => "HC",
            Self::Light => "LC",
        }
    }

    fn from_prefix(value: &str) -> Option<Self> {
        match value {
            "HC" => Some(Self::Heavy),
            "LC" => Some(Self::Light),
            _ => None,
        }
    }
}

/// Compact, index-resolved structural identity for one complete V(D)J event.
///
/// Segment identities are stored as 12-bit indices into the persisted VDJ
/// reference. Common measurements 0..14 use one hexadecimal digit; larger u16
/// values use `F` followed by four hexadecimal digits. The prefix describes the
/// recombination geometry: HC carries V-D-J, LC carries V-J. The actual IG/TR
/// locus is recovered from the referenced segments and validated on decode.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PackedRecombinationId {
    role: ReceptorRole,
    payload: String,
}

impl PackedRecombinationId {
    const VERSION: u8 = 1;

    pub fn from_call(call: &RearrangementCall) -> Option<Self> {
        let v = call.v.as_ref()?;
        let j = call.j.as_ref()?;
        let junction = call.junction.as_ref()?;
        if v.segment_index > 0x0fff || j.segment_index > 0x0fff {
            return None;
        }
        let role = if call.chain.has_d() {
            ReceptorRole::Heavy
        } else {
            ReceptorRole::Light
        };
        let mut payload = format!("{:X}{:03X}", Self::VERSION, v.segment_index);
        if role == ReceptorRole::Heavy {
            let d = call.d.as_ref()?;
            if d.segment_index > 0x0fff {
                return None;
            }
            payload.push_str(&format!("{:03X}", d.segment_index));
        }
        payload.push_str(&format!("{:03X}", j.segment_index));
        push_measurement(&mut payload, junction.v_del_3);
        push_measurement(&mut payload, junction.p_v3_len());
        push_measurement(&mut payload, junction.n1_len());
        if role == ReceptorRole::Heavy {
            push_measurement(&mut payload, junction.p_d5_len());
            push_measurement(&mut payload, junction.d_del_5?);
            push_measurement(&mut payload, junction.d_retained_len?);
            push_measurement(&mut payload, junction.d_del_3?);
            push_measurement(&mut payload, junction.p_d3_len());
            push_measurement(&mut payload, junction.n2_len());
        }
        push_measurement(&mut payload, junction.p_j5_len());
        push_measurement(&mut payload, junction.j_del_5);
        payload.push(if junction.pn_alternative { '1' } else { '0' });
        Some(Self { role, payload })
    }

    pub fn role(&self) -> ReceptorRole {
        self.role
    }

    pub fn payload(&self) -> &str {
        &self.payload
    }

    pub fn decode(&self, reference: &crate::reference::VdjReference) -> Result<RecombinationMeasurements, String> {
        let mut cursor = HexCursor::new(&self.payload);
        let version = cursor.read_nibble()?;
        if version != Self::VERSION {
            return Err(format!("unsupported {} encoding version {version}", self.role.prefix()));
        }
        let v_index = cursor.read_segment_index()?;
        let d_index = if self.role == ReceptorRole::Heavy {
            Some(cursor.read_segment_index()?)
        } else {
            None
        };
        let j_index = cursor.read_segment_index()?;

        let v = reference.segments.get(v_index).ok_or_else(|| format!("V segment index {v_index} is outside this VDJ index"))?;
        let j = reference.segments.get(j_index).ok_or_else(|| format!("J segment index {j_index} is outside this VDJ index"))?;
        if v.kind != crate::types::SegmentKind::V {
            return Err(format!("segment {v_index} ({}) is not a V segment", v.name));
        }
        if j.kind != crate::types::SegmentKind::J {
            return Err(format!("segment {j_index} ({}) is not a J segment", j.name));
        }
        if v.chain != j.chain {
            return Err(format!("V {} and J {} resolve to different chains", v.name, j.name));
        }
        if v.chain.has_d() != (self.role == ReceptorRole::Heavy) {
            return Err(format!("{} prefix is incompatible with resolved chain {}", self.role.prefix(), v.chain));
        }
        let d = if let Some(d_index) = d_index {
            let d = reference.segments.get(d_index).ok_or_else(|| format!("D segment index {d_index} is outside this VDJ index"))?;
            if d.kind != crate::types::SegmentKind::D || d.chain != v.chain {
                return Err(format!("segment {d_index} ({}) is not a matching {} D segment", d.name, v.chain));
            }
            Some(d)
        } else {
            None
        };

        let v_del_3 = Some(cursor.read_measurement()?);
        let p_v3_len = Some(cursor.read_measurement()?);
        let n1_len = Some(cursor.read_measurement()?);
        let (p_d5_len, d_del_5, d_retained_len, d_del_3, p_d3_len, n2_len) = if self.role == ReceptorRole::Heavy {
            (
                Some(cursor.read_measurement()?),
                Some(cursor.read_measurement()?),
                Some(cursor.read_measurement()?),
                Some(cursor.read_measurement()?),
                Some(cursor.read_measurement()?),
                Some(cursor.read_measurement()?),
            )
        } else {
            (None, None, None, None, None, None)
        };
        let p_j5_len = Some(cursor.read_measurement()?);
        let j_del_5 = Some(cursor.read_measurement()?);
        let pn_alternative = match cursor.read_nibble()? {
            0 => false,
            1 => true,
            other => return Err(format!("invalid P/N alternative flag {other:X}")),
        };
        if !cursor.is_finished() {
            return Err("trailing hexadecimal data in recombination id".into());
        }

        Ok(RecombinationMeasurements {
            chain: v.chain,
            v: v.name.clone(),
            v_del_3,
            p_v3_len,
            n1_len,
            p_d5_len,
            d: d.map(|x| x.name.clone()),
            d_del_5,
            d_retained_len,
            d_del_3,
            p_d3_len,
            n2_len,
            p_j5_len,
            j_del_5,
            j: j.name.clone(),
            pn_alternative,
        })
    }
}

impl fmt::Display for PackedRecombinationId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.role.prefix(), self.payload)
    }
}

impl FromStr for PackedRecombinationId {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let (prefix, payload) = text.split_once(':').ok_or_else(|| "packed recombination id requires HC: or LC: prefix".to_string())?;
        let role = ReceptorRole::from_prefix(prefix).ok_or_else(|| format!("unknown recombination id prefix: {prefix}"))?;
        if payload.is_empty() || !payload.bytes().all(|x| x.is_ascii_hexdigit()) {
            return Err("packed recombination id payload must be hexadecimal".into());
        }
        Ok(Self { role, payload: payload.to_ascii_uppercase() })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LightChainStatus {
    Present,
    PossiblyRecombining,
    NotDetected,
}

impl fmt::Display for LightChainStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Present => "present",
            Self::PossiblyRecombining => "possibly_recombining",
            Self::NotDetected => "not_detected",
        })
    }
}

impl CellVdjSummary {
    pub fn strongest_heavy_id(&self) -> Option<PackedRecombinationId> {
        self.rearrangements
            .iter()
            .filter(|call| call.chain.has_d())
            .max_by(|a, b| compare_calls(a, b))
            .and_then(PackedRecombinationId::from_call)
    }

    pub fn strongest_light_id(&self) -> Option<PackedRecombinationId> {
        self.rearrangements
            .iter()
            .filter(|call| !call.chain.has_d())
            .max_by(|a, b| compare_calls(a, b))
            .and_then(PackedRecombinationId::from_call)
    }

    pub fn light_chain_status(&self) -> LightChainStatus {
        if self.strongest_light_id().is_some() {
            LightChainStatus::Present
        } else if self.recombination_activity.rag_pair_detected {
            LightChainStatus::PossiblyRecombining
        } else {
            LightChainStatus::NotDetected
        }
    }
}

fn push_measurement(out: &mut String, value: u16) {
    if value <= 14 {
        out.push(char::from_digit(u32::from(value), 16).unwrap().to_ascii_uppercase());
    } else {
        out.push('F');
        out.push_str(&format!("{value:04X}"));
    }
}

struct HexCursor<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> HexCursor<'a> {
    fn new(text: &'a str) -> Self {
        Self { bytes: text.as_bytes(), pos: 0 }
    }

    fn read_nibble(&mut self) -> Result<u8, String> {
        let byte = *self.bytes.get(self.pos).ok_or_else(|| "truncated recombination id".to_string())?;
        self.pos += 1;
        hex_nibble(byte).ok_or_else(|| format!("invalid hex digit {}", byte as char))
    }

    fn read_segment_index(&mut self) -> Result<usize, String> {
        let mut value = 0usize;
        for _ in 0..3 {
            value = (value << 4) | usize::from(self.read_nibble()?);
        }
        Ok(value)
    }

    fn read_measurement(&mut self) -> Result<u16, String> {
        let first = self.read_nibble()?;
        if first < 15 {
            return Ok(u16::from(first));
        }
        let mut value = 0u16;
        for _ in 0..4 {
            value = (value << 4) | u16::from(self.read_nibble()?);
        }
        Ok(value)
    }

    fn is_finished(&self) -> bool {
        self.pos == self.bytes.len()
    }
}

fn hex_nibble(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}

fn compare_calls(a: &RearrangementCall, b: &RearrangementCall) -> std::cmp::Ordering {
    let score = |call: &RearrangementCall| {
        call.v.as_ref().map_or(0, |x| x.local_alignment_score)
            + call.j.as_ref().map_or(0, |x| x.local_alignment_score)
            + call.d.as_ref().map_or(0, |x| x.local_alignment_score)
    };
    a.total_supporting_umis
        .cmp(&b.total_supporting_umis)
        .then_with(|| score(a).cmp(&score(b)))
        .then_with(|| a.chain.cmp(&b.chain))
        .then_with(|| a.notation.cmp(&b.notation))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::junction::JunctionMeasurement;
    use crate::posterior::{GermlineSegmentSupport, RearrangementCall, RecombinationStage};
    use crate::reference::VdjReference;
    use crate::types::{SegmentKind, VdjSegment};

    fn support(id: &str, kind: SegmentKind, score: i32, segment_index: usize) -> GermlineSegmentSupport {
        GermlineSegmentSupport {
            segment_index,
            id: id.into(),
            kind,
            local_alignment_score: score,
            supporting_umis: 1,
            supporting_reads: 1,
            locus_fraction: 0.0,
            distance_to_recombination_center: 0,
        }
    }

    fn call(chain: Chain, v: &str, d: Option<&str>, j: &str, umis: usize, score: i32) -> RearrangementCall {
        RearrangementCall {
            chain,
            stage: if chain.has_d() { RecombinationStage::Vdj } else { RecombinationStage::Vj },
            v: Some(support(v, SegmentKind::V, score, 0)),
            d: d.map(|id| support(id, SegmentKind::D, score, 1)),
            d_inferred_from_vj_junction: false,
            d_hypothesis_margin: None,
            j: Some(support(j, SegmentKind::J, score, if d.is_some() { 2 } else { 1 })),
            c: None,
            total_supporting_umis: umis,
            supporting_reads: Vec::new(),
            junction: None,
            notation: format!("{chain}:{v}-{j}"),
        }
    }

    fn junction(vdj: bool) -> JunctionMeasurement {
        JunctionMeasurement {
            v_del_3: 0,
            p_v3: b"AAA".to_vec(),
            n1: b"C".to_vec(),
            p_d5: Vec::new(),
            d_del_5: vdj.then_some(0),
            d_retained_len: vdj.then_some(23),
            d_del_3: vdj.then_some(0),
            p_d3: if vdj { b"G".to_vec() } else { Vec::new() },
            n2: if vdj { b"T".to_vec() } else { Vec::new() },
            p_j5: Vec::new(),
            j_del_5: 4,
            pn_alternative: vdj,
            observed_v: b"ACGT".to_vec(),
            observed_d: if vdj { b"GG".to_vec() } else { Vec::new() },
            observed_j: b"TTAA".to_vec(),
            naive_v: b"ACGT".to_vec(),
            naive_d: if vdj { b"GG".to_vec() } else { Vec::new() },
            naive_j: b"TTAA".to_vec(),
            observed_sequence: Vec::new(),
            inferred_naive_sequence: Vec::new(),
        }
    }

    #[test]
    fn packed_hc_is_short_and_decodes_through_reference() {
        let segment = |name: &str, kind: SegmentKind| VdjSegment {
            name: name.into(), transcript_id: format!("{name}.tx"), gene_id: format!("{name}.gene"),
            chain: Chain::Igh, kind, chr: "chr12".into(), start: 0, end: 10,
            strand_minus: false, locus_rank: 0, locus_fraction: 0.0,
            distance_to_recombination_center: 0, sequence: b"ACGTACGTAC".to_vec(),
        };
        let reference = VdjReference { segments: vec![
            segment("Ighv1-64", SegmentKind::V),
            segment("Ighd1-1", SegmentKind::D),
            segment("Ighj3", SegmentKind::J),
        ]};
        let mut call = call(Chain::Igh, "Ighv1-64", Some("Ighd1-1"), "Ighj3", 2, 100);
        call.junction = Some(junction(true));
        let id = PackedRecombinationId::from_call(&call).unwrap();
        assert!(id.to_string().starts_with("HC:"));
        assert!(id.to_string().len() < 40);
        assert!(id.payload().contains("F0017"));
        let decoded = id.decode(&reference).unwrap();
        assert_eq!(decoded.chain, Chain::Igh);
        assert_eq!(decoded.v, "Ighv1-64");
        assert_eq!(decoded.d.as_deref(), Some("Ighd1-1"));
        assert_eq!(decoded.j, "Ighj3");
        assert_eq!(decoded.d_retained_len, Some(23));
        assert!(decoded.pn_alternative);
    }

    #[test]
    fn packed_lc_uses_vj_geometry() {
        let mut call = call(Chain::Igk, "Igkv6-15", None, "Igkj2", 9, 100);
        call.v.as_mut().unwrap().segment_index = 12;
        call.j.as_mut().unwrap().segment_index = 34;
        call.junction = Some(junction(false));
        let id = PackedRecombinationId::from_call(&call).unwrap();
        assert!(id.to_string().starts_with("LC:"));
        assert_eq!(id.to_string().len(), 16);
    }
}
