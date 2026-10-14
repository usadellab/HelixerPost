use std::fmt::{Display, Formatter};
use std::io::{BufWriter, Write};

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum GffFeature {
    Gene,
    MRNA,
    Exon,
    FivePrimeUTR,
    CDS,
    ThreePrimeUTR,
}

impl GffFeature {
    pub fn as_str(self) -> &'static str {
        match self {
            GffFeature::Gene => "gene",
            GffFeature::MRNA => "mRNA",
            GffFeature::Exon => "exon",
            GffFeature::FivePrimeUTR => "five_prime_UTR",
            GffFeature::CDS => "CDS",
            GffFeature::ThreePrimeUTR => "three_prime_UTR",
        }
    }
}

impl Display for GffFeature {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum GffStrand {
    Forward,
    Reverse,
}

impl GffStrand {
    pub fn other(self) -> GffStrand {
        match self {
            GffStrand::Forward => GffStrand::Reverse,
            GffStrand::Reverse => GffStrand::Forward,
        }
    }

    pub fn from_rev(rev: bool) -> GffStrand {
        if rev {
            GffStrand::Reverse
        } else {
            GffStrand::Forward
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            GffStrand::Forward => "+",
            GffStrand::Reverse => "-",
        }
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum GffPhase {
    Zero = 0,
    One = 1,
    Two = 2,
}

impl GffPhase {
    pub fn as_str(self) -> &'static str {
        match self {
            GffPhase::Zero => "0",
            GffPhase::One => "1",
            GffPhase::Two => "2",
        }
    }
}

impl From<u64> for GffPhase {
    fn from(offset: u64) -> Self {
        match offset % 3 {
            0 => GffPhase::Zero,
            1 => GffPhase::Two, // Phase means 'bases left to read', so counts down - Because biology
            2 => GffPhase::One,
            _ => panic!("Math is broken"),
        }
    }
}

/*

1	sequence	The name of the sequence where the feature is located.
2	source	Keyword identifying the source of the feature, like a program (e.g. Augustus or RepeatMasker) or an organization (like TAIR).
3	feature	The feature type name, like "gene" or "exon". In a well structured GFF file, all the children features always follow their parents in a single block (so all exons of a transcript are put after their parent "transcript" feature line and before any other parent transcript line). In GFF3, all features and their relationships should be compatible with the standards released by the Sequence Ontology Project.
4	start	Genomic start of the feature, with a 1-base offset. This is in contrast with other 0-offset half-open sequence formats, like BED.
5	end	Genomic end of the feature, with a 1-base offset. This is the same end coordinate as it is in 0-offset half-open sequence formats, like BED.[citation needed]
6	score	Numeric value that generally indicates the confidence of the source in the annotated feature. A value of "." (a dot) is used to define a null value.
7	strand	Single character that indicates the strand of the feature; it can assume the values of "+" (positive, or 5'->3'), "-", (negative, or 3'->5'), "." (undetermined).
8	phase	phase of CDS features; it can be either one of 0, 1, 2 (for CDS features) or "." (for everything else). See the section below for a detailed explanation.
9	attributes	All the other information pertaining to this feature. The format, structure and content of this field is the one which varies the most between the three competing file formats.

 */

pub struct GffRecord {
    sequence: String,
    source: String,
    feature: GffFeature,
    start: u64,
    end: u64,
    score: Option<f32>,
    strand: Option<GffStrand>,
    phase: Option<GffPhase>,
    attributes: String,
}

impl GffRecord {
    pub fn new(
        sequence: String,
        source: String,
        feature: GffFeature,
        start: u64,
        end: u64,
        score: Option<f32>,
        strand: Option<GffStrand>,
        phase: Option<GffPhase>,
        attributes: String,
    ) -> GffRecord {
        GffRecord {
            sequence,
            source,
            feature,
            start,
            end,
            score,
            strand,
            phase,
            attributes,
        }
    }

    pub fn swap_strand(&mut self, len: u64) {
        let start = self.start;
        let end = self.end;

        self.end = 1 + len - start;
        self.start = 1 + len - end;

        if let Some(strand) = self.strand {
            self.strand = Some(strand.other())
        };
    }

    pub fn get_sequence(&self) -> &String {
        &self.sequence
    }

    pub fn get_source(&self) -> &String {
        &self.source
    }

    pub fn get_start(&self) -> u64 {
        self.start
    }

    pub fn get_end(&self) -> u64 {
        self.end
    }

    pub fn get_score(&self) -> Option<f32> {
        self.score
    }

    pub fn get_strand(&self) -> Option<GffStrand> {
        self.strand
    }

    pub fn get_phase(&self) -> Option<GffPhase> {
        self.phase
    }

    pub fn get_attributes(&self) -> &String {
        &self.attributes
    }
}

pub struct GffWriter<W: Write> {
    writer: BufWriter<W>,
}

impl<W: Write> GffWriter<W> {
    pub fn new(writer: BufWriter<W>) -> GffWriter<W> {
        GffWriter { writer }
    }

    /// Writes header for the top of the GFF3 file
    /// and optionally includes species name and Helixer model info
    ///
    /// # Arguments
    ///
    /// * `species` - An Option containing the species name
    /// * `helixer_model_md5sum` - An Option containing the md5checksum
    /// and file path information from Helixer (as found in h5.attrs\["model_md5sum"\] in the predictions output.)
    pub fn write_global_header(
        &mut self,
        species: Option<&str>,
        helixer_model_md5sum: Option<String>,
    ) -> std::io::Result<()> {
        const GFF_VERSION: &'static str = "3.2.1";
        write!(self.writer, "##gff-version {}\n", GFF_VERSION)?;
        if let Some(species) = species {
            write!(self.writer, "##species {}\n", species)?;
        }
        if let Some(helixer_model_md5sum) = helixer_model_md5sum {
            write!(self.writer, "# {}\n", helixer_model_md5sum)?;
        }
        Ok(())
    }

    /// Writes header for each region in the GFF3 file indicating the name and boundaries of the region. For Helixer output, regions always start at 1.
    ///
    /// # Arguments
    ///
    /// * `sequence_name` - name of the sequence aka region
    /// * `sequence length` - length of the sequence aka region, since sequences start at 1 this will also be the end of the region.
    pub fn write_region_header(
        &mut self,
        sequence_name: &str,
        sequence_length: u64,
    ) -> std::io::Result<()> {
        let region_info = format!("{} {} {}", sequence_name, 1, sequence_length);
        write!(self.writer, "##sequence-region {}\n", region_info)
    }

    pub fn write_record(&mut self, rec: &GffRecord) -> std::io::Result<()> {
        let score = rec.score.map_or(".".to_owned(), |v| format!("{}", v));
        let strand = rec.strand.map_or(".", |v| v.as_str());
        let phase = rec.phase.map_or(".", |v| v.as_str());

        write!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            rec.sequence,
            rec.source,
            rec.feature,
            rec.start,
            rec.end,
            score,
            strand,
            phase,
            rec.attributes
        )
    }

    pub fn write_records(&mut self, recs: &[GffRecord]) -> std::io::Result<()> {
        for rec in recs {
            self.write_record(rec)?
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(start: u64, end: u64, strand: GffStrand) -> GffRecord {
        GffRecord::new(
            "seq".to_owned(),
            "src".to_owned(),
            GffFeature::Exon,
            start,
            end,
            None,
            Some(strand),
            None,
            String::new(),
        )
    }

    /// The `From<u64> for GffPhase` mapping counts backwards because the GFF
    /// "phase" column is bases-left-to-read, not bases-consumed. This test
    /// pins that semantics: any future "simplification" that swaps One/Two
    /// will fail loudly.
    #[test]
    fn gff_phase_counts_down_because_biology() {
        assert_eq!(GffPhase::from(0_u64), GffPhase::Zero);
        assert_eq!(GffPhase::from(1_u64), GffPhase::Two);
        assert_eq!(GffPhase::from(2_u64), GffPhase::One);
        assert_eq!(GffPhase::from(3_u64), GffPhase::Zero);
        assert_eq!(GffPhase::from(7_u64), GffPhase::Two);
    }

    /// swap_strand is invoked once per gene on the reverse strand. A wrong
    /// formula would silently shift every rev-strand annotation.
    #[test]
    fn swap_strand_known_position() {
        // Sequence length 1000, feature at 100..=200 on the forward strand.
        // Reverse-complement coordinates: start = 1+1000-200 = 801,
        // end = 1+1000-100 = 901.
        let mut r = rec(100, 200, GffStrand::Forward);
        r.swap_strand(1000);
        assert_eq!(r.get_start(), 801);
        assert_eq!(r.get_end(), 901);
        assert_eq!(r.get_strand(), Some(GffStrand::Reverse));
    }

    /// A double swap_strand must return the original coordinates and strand.
    #[test]
    fn swap_strand_is_an_involution() {
        let mut r = rec(120, 480, GffStrand::Forward);
        r.swap_strand(2000);
        r.swap_strand(2000);
        assert_eq!(r.get_start(), 120);
        assert_eq!(r.get_end(), 480);
        assert_eq!(r.get_strand(), Some(GffStrand::Forward));
    }

    /// A feature spanning the entire sequence must map to itself under
    /// reverse-complement.
    #[test]
    fn swap_strand_at_full_sequence_boundary() {
        let mut r = rec(1, 100, GffStrand::Forward);
        r.swap_strand(100);
        assert_eq!(r.get_start(), 1);
        assert_eq!(r.get_end(), 100);
        assert_eq!(r.get_strand(), Some(GffStrand::Reverse));
    }
}
