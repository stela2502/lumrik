
use std::fmt;

#[derive(Debug, Clone )]
pub struct BamFlag {
    flag: u16
}

impl Default for BamFlag {
    fn default() -> Self {
        Self {
            flag: 0,
        }
    }
}

// Implementing Display trait for BamFlag
impl fmt::Display for BamFlag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.flag )
    }
}

impl BamFlag {

        /// Create a new BamFlags object from a numeric value
    pub fn new(flag: u16) -> Self {
        Self { flag }
    }

    /// Check if the flag indicates paired reads
    pub fn is_paired(&self) -> bool {
        self.flag & 0x1 != 0
    }

    /// Check if the flag indicates proper pair
    pub fn is_proper_pair(&self) -> bool {
        self.flag & 0x2 != 0
    }

    /// Check if the flag indicates unmapped read
    pub fn is_unmapped(&self) -> bool {
        self.flag & 0x4 != 0
    }

    /// Check if the flag indicates mate unmapped
    pub fn is_mate_unmapped(&self) -> bool {
        self.flag & 0x8 != 0
    }

    /// Check if the flag indicates reverse strand
    pub fn is_reverse_strand(&self) -> bool {
        self.flag & 0x10 != 0
    }

    /// Check if the flag indicates mate reverse strand
    pub fn is_mate_reverse_strand(&self) -> bool {
        self.flag & 0x20 != 0
    }

    /// Check if the flag indicates read 1
    pub fn is_read1(&self) -> bool {
        self.flag & 0x40 != 0
    }

    /// Check if the flag indicates read 2
    pub fn is_read2(&self) -> bool {
        self.flag & 0x80 != 0
    }

    /// Check if the flag indicates secondary alignment
    pub fn is_secondary(&self) -> bool {
        self.flag & 0x100 != 0
    }

    /// Check if the flag indicates a failed quality check
    pub fn is_qc_fail(&self) -> bool {
        self.flag & 0x200 != 0
    }

    /// Check if the flag indicates a duplicate read
    pub fn is_duplicate(&self) -> bool {
        self.flag & 0x400 != 0
    }  


    /// Convert the boolean flags into the corresponding SAM flag integer
    pub fn to_sam(&self) -> u16 {
        self.flag
    }

    /// simple way to check the deffernet values like flag.is("paired")
    pub fn is(&self, tag: &str) -> bool {
        match tag {
            "paired" => self.is_paired(),
            "proper_pair" => self.is_proper_pair(),
            "unmapped" => self.is_unmapped(),
            "mate_unmapped" => self.is_mate_unmapped(),
            "reverse_strand" => self.is_reverse_strand(),
            "mate_reverse_strand" => self.is_mate_reverse_strand(),
            "read1" => self.is_read1(),
            "read2" => self.is_read2(),
            "secondary" => self.is_secondary(),
            "qc_fail" => self.is_qc_fail(),
            "duplicate" => self.is_duplicate(),
            _ => false, // Return false for unknown tags
        }
    }


        // Setter methods to modify specific flags
    pub fn set_paired(&mut self, paired: bool) {
        if paired {
            self.flag |= 0x1; // Set the paired flag (bit 0)
        } else {
            self.flag &= !0x1; // Clear the paired flag (bit 0)
        }
    }

    pub fn set_proper_pair(&mut self, proper_pair: bool) {
        if proper_pair {
            self.flag |= 0x2; // Set the proper pair flag (bit 1)
        } else {
            self.flag &= !0x2; // Clear the proper pair flag (bit 1)
        }
    }

    pub fn set_unmapped(&mut self, unmapped: bool) {
        if unmapped {
            self.flag |= 0x4; // Set the unmapped flag (bit 2)
        } else {
            self.flag &= !0x4; // Clear the unmapped flag (bit 2)
        }
    }

    pub fn set_mate_unmapped(&mut self, mate_unmapped: bool) {
        if mate_unmapped {
            self.flag |= 0x8; // Set the mate unmapped flag (bit 3)
        } else {
            self.flag &= !0x8; // Clear the mate unmapped flag (bit 3)
        }
    }

    pub fn set_reverse_strand(&mut self, reverse_strand: bool) {
        if reverse_strand {
            self.flag |= 0x10; // Set the reverse strand flag (bit 4)
        } else {
            self.flag &= !0x10; // Clear the reverse strand flag (bit 4)
        }
    }

    pub fn set_mate_reverse_strand(&mut self, mate_reverse_strand: bool) {
        if mate_reverse_strand {
            self.flag |= 0x20; // Set the mate reverse strand flag (bit 5)
        } else {
            self.flag &= !0x20; // Clear the mate reverse strand flag (bit 5)
        }
    }

    pub fn set_read1(&mut self, read1: bool) {
        if read1 {
            self.flag |= 0x40; // Set the read1 flag (bit 6)
        } else {
            self.flag &= !0x40; // Clear the read1 flag (bit 6)
        }
    }

    pub fn set_read2(&mut self, read2: bool) {
        if read2 {
            self.flag |= 0x80; // Set the read2 flag (bit 7)
        } else {
            self.flag &= !0x80; // Clear the read2 flag (bit 7)
        }
    }

    pub fn set_secondary(&mut self, secondary: bool) {
        if secondary {
            self.flag |= 0x100; // Set the secondary alignment flag (bit 8)
        } else {
            self.flag &= !0x100; // Clear the secondary alignment flag (bit 8)
        }
    }

    pub fn set_qc_fail(&mut self, qc_fail: bool) {
        if qc_fail {
            self.flag |= 0x200; // Set the QC fail flag (bit 9)
        } else {
            self.flag &= !0x200; // Clear the QC fail flag (bit 9)
        }
    }

    pub fn set_duplicate(&mut self, duplicate: bool) {
        if duplicate {
            self.flag |= 0x400; // Set the duplicate flag (bit 10)
        } else {
            self.flag &= !0x400; // Clear the duplicate flag (bit 10)
        }
    }

}

