use crate::error::{PrimerError, PrimerResult};
use std::fmt;
use std::ops::Range;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Forward,
    ReverseComplement,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimerSlice<'a> {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    source_seq: Option<&'a [u8]>,
    source_qual: Option<&'a [u8]>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimerSegment {
    pub name: String,
    pub ranges: Vec<Range<usize>>,
}

impl fmt::Display for PrimerSegment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: ", self.name)?;

        for (idx, range) in self.ranges.iter().enumerate() {
            if idx > 0 {
                write!(f, ", ")?;
            }

            write!(f, "{}..{}", range.start, range.end)?;
        }

        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimerMatch {
    pub chemistry_name: String,
    pub orientation: Orientation,
    pub primer_start: usize,
    pub primer_end: usize,
    pub insert_start: usize,
    pub insert_end: usize,
    pub bd_cell_id: Option<u64>,
    pub cell_seq: Option<Vec<u8>>,
    pub segments: Vec<PrimerSegment>,
}

impl fmt::Display for PrimerMatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "{} {:?} primer={}..{} insert={}..{}",
            self.chemistry_name,
            self.orientation,
            self.primer_start,
            self.primer_end,
            self.insert_start,
            self.insert_end,
        )?;

        if let Some(id) = self.bd_cell_id {
            writeln!(f, "  bd_cell_id={id}")?;
        }

        if let Some(seq) = &self.cell_seq {
            writeln!(f, "  cell_seq={}", String::from_utf8_lossy(seq))?;
        }

        for segment in &self.segments {
            writeln!(f, "  {segment}")?;
        }

        Ok(())
    }
}

impl<'a> PrimerSlice<'a> {
    pub fn new(seq: Vec<u8>, qual: Vec<u8>) -> Self {
        Self {
            seq,
            qual,
            source_seq: None,
            source_qual: None,
        }
    }

    pub fn from_range(seq: &'a [u8], qual: &'a [u8], range: Range<usize>) -> PrimerResult<Self> {
        Self::check_inputs(seq, qual, std::slice::from_ref(&range))?;
        Ok(Self {
            seq: seq[range.clone()].to_vec(),
            qual: qual[range].to_vec(),
            source_seq: None,
            source_qual: None,
        })
    }

    pub fn from_ranges(
        seq: &'a [u8],
        qual: &'a [u8],
        ranges: &[Range<usize>],
    ) -> PrimerResult<Self> {
        Self::check_inputs(seq, qual, ranges)?;
        let len = ranges.iter().map(|range| range.end - range.start).sum();
        let mut out_seq = Vec::with_capacity(len);
        let mut out_qual = Vec::with_capacity(len);
        for range in ranges {
            out_seq.extend_from_slice(&seq[range.clone()]);
            out_qual.extend_from_slice(&qual[range.clone()]);
        }
        Ok(Self {
            seq: out_seq,
            qual: out_qual,
            source_seq: None,
            source_qual: None,
        })
    }

    pub fn check_inputs(seq: &[u8], qual: &[u8], ranges: &[Range<usize>]) -> PrimerResult<()> {
        if seq.len() != qual.len() {
            return Err(PrimerError::invalid_coordinates(
                "sequence and quality have different lengths",
            ));
        }
        for range in ranges {
            if range.start > range.end || range.end > seq.len() {
                return Err(PrimerError::invalid_coordinates(format!(
                    "range {}..{} outside sequence length {}",
                    range.start,
                    range.end,
                    seq.len()
                )));
            }
        }
        Ok(())
    }
}

impl PrimerSegment {
    pub fn new(name: impl Into<String>, ranges: Vec<Range<usize>>) -> Self {
        Self {
            name: name.into(),
            ranges,
        }
    }
}

impl PrimerMatch {
    pub fn new(chemistry_name: String, orientation: Orientation) -> Self {
        Self {
            chemistry_name,
            orientation,
            primer_start: 0,
            primer_end: 0,
            insert_start: 0,
            insert_end: 0,
            bd_cell_id: None,
            cell_seq: None,
            segments: Vec::new(),
        }
    }

    pub fn add_cell_seq(&mut self, cell_seq: impl AsRef<[u8]>) {
        self.cell_seq = Some(cell_seq.as_ref().to_vec());
    }

    pub fn add_segment(&mut self, name: impl Into<String>, range: Range<usize>) {
        self.segments.push(PrimerSegment::new(name, vec![range]));
    }

    pub fn add_segment_ranges(&mut self, name: &str, ranges: Vec<Range<usize>>) {
        self.segments.push(PrimerSegment::new(name, ranges));
    }

    pub fn get_cell<'a>(&self, seq: &'a [u8], qual: &'a [u8]) -> PrimerResult<PrimerSlice<'a>> {
        self.get_segment(seq, qual, "CELL")
            .or_else(|_| self.get_segment(seq, qual, "BD_CELL"))
    }

    pub fn get_umi<'a>(&self, seq: &'a [u8], qual: &'a [u8]) -> PrimerResult<PrimerSlice<'a>> {
        self.get_segment(seq, qual, "UMI")
    }

    pub fn get_insert<'a>(&self, seq: &'a [u8], qual: &'a [u8]) -> PrimerResult<PrimerSlice<'a>> {
        PrimerSlice::from_range(seq, qual, self.insert_start..self.insert_end)
    }

    pub fn get_segment<'a>(
        &self,
        seq: &'a [u8],
        qual: &'a [u8],
        name: &str,
    ) -> PrimerResult<PrimerSlice<'a>> {
        let Some(segment) = self.segments.iter().find(|segment| segment.name == name) else {
            return Err(PrimerError::invalid_coordinates(format!(
                "segment '{name}' is not present"
            )));
        };
        PrimerSlice::from_ranges(seq, qual, &segment.ranges)
    }

    pub fn shift_by(&mut self, offset: usize) {
        self.primer_start += offset;
        self.primer_end += offset;
        self.insert_start += offset;
        self.insert_end += offset;
        for segment in &mut self.segments {
            for range in &mut segment.ranges {
                range.start += offset;
                range.end += offset;
            }
        }
    }

    pub fn remap_reverse_coordinates(&mut self, len: usize) {
        self.primer_start = len - self.primer_start;
        self.primer_end = len - self.primer_end;
        std::mem::swap(&mut self.primer_start, &mut self.primer_end);
        self.insert_start = len - self.insert_start;
        self.insert_end = len - self.insert_end;
        std::mem::swap(&mut self.insert_start, &mut self.insert_end);

        for segment in &mut self.segments {
            for range in &mut segment.ranges {
                let old_start = range.start;
                let old_end = range.end;
                range.start = len - old_end;
                range.end = len - old_start;
            }
            segment.ranges.reverse();
        }
    }
}

#[derive(Debug, Clone)]
pub struct PrimerAttempt {
    pub offset: usize,
    pub orientation: Orientation,
    pub ok: bool,
    pub reason: String,
    pub cell_seq: Option<String>,
    pub segments: Vec<PrimerSegmentAttempt>,
}

#[derive(Debug, Clone)]
pub struct PrimerSegmentAttempt {
    pub name: String,
    pub range: std::ops::Range<usize>,
    pub dna: String,
    pub ok: bool,
    pub reason: String,
}
