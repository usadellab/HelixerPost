/*

Chr1    phytozomev10    gene    91376   95651   .       +       .       ID=AT1G01220.TAIR10;Name=AT1G01220
Chr1    phytozomev10    mRNA    91376   95651   .       +       .       ID=AT1G01220.1.TAIR10;Name=AT1G01220.1;pacid=19652890;longest=1;Parent=AT1G01220.TAIR10
Chr1    phytozomev10    exon    91376   91633   .       +       .       ID=AT1G01220.1.TAIR10.exon.1;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    five_prime_UTR  91376   91633   .       +       .       ID=AT1G01220.1.TAIR10.five_prime_UTR.1;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    exon    91743   92070   .       +       .       ID=AT1G01220.1.TAIR10.exon.2;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    five_prime_UTR  91743   91749   .       +       .       ID=AT1G01220.1.TAIR10.five_prime_UTR.2;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    CDS     91750   92070   .       +       0       ID=AT1G01220.1.TAIR10.CDS.1;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    exon    92270   92501   .       +       .       ID=AT1G01220.1.TAIR10.exon.3;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    CDS     92270   92501   .       +       0       ID=AT1G01220.1.TAIR10.CDS.2;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    exon    92569   92933   .       +       .       ID=AT1G01220.1.TAIR10.exon.4;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    CDS     92569   92933   .       +       2       ID=AT1G01220.1.TAIR10.CDS.3;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    exon    93045   93171   .       +       .       ID=AT1G01220.1.TAIR10.exon.5;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    CDS     93045   93171   .       +       0       ID=AT1G01220.1.TAIR10.CDS.4;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    exon    93271   94281   .       +       .       ID=AT1G01220.1.TAIR10.exon.6;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    CDS     93271   94281   .       +       2       ID=AT1G01220.1.TAIR10.CDS.5;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    exon    94357   95075   .       +       .       ID=AT1G01220.1.TAIR10.exon.7;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    CDS     94357   95075   .       +       2       ID=AT1G01220.1.TAIR10.CDS.6;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    exon    95160   95651   .       +       .       ID=AT1G01220.1.TAIR10.exon.8;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    CDS     95160   95552   .       +       0       ID=AT1G01220.1.TAIR10.CDS.7;Parent=AT1G01220.1.TAIR10;pacid=19652890
Chr1    phytozomev10    three_prime_UTR 95553   95651   .       +       .       ID=AT1G01220.1.TAIR10.three_prime_UTR.1;Parent=AT1G01220.1.TAIR10;pacid=19652890


sequence: String, source: String, feature: GffFeature, start: u64, end: u64, score: Option<f32>,
           strand: Option<GffStrand>, phase: Option<GffPhase>, attributes: String

 */

use crate::analysis::hmm::{HmmAnnotationLabel, HmmStateRegion};
use crate::gff::{GffFeature, GffPhase, GffRecord, GffStrand};

// generate gene, mRNA and exon records, based on UTR5/CDS/UTR3
fn generate_gff_aggregate_records(
    recs: Vec<GffRecord>,
    sequence: &str,
    source: &str,
    strand: Option<GffStrand>,
    gene_name: &str,
) -> Vec<GffRecord> {
    if recs.len() == 0 {
        return Vec::new();
    }

    let mut maybe_transcript_start: Option<u64> = None;
    let mut maybe_transcript_end: Option<u64> = None;

    let mut maybe_exon_start: Option<u64> = None;
    let mut maybe_exon_end: Option<u64> = None;

    let mut exon_ranges = Vec::new();

    for rec in recs.iter() {
        if maybe_transcript_start == None {
            maybe_transcript_start = Some(rec.get_start());
        }
        maybe_transcript_end = Some(rec.get_end());

        if let (Some(exon_start), Some(exon_end)) = (maybe_exon_start, maybe_exon_end) {
            if exon_end + 1 < rec.get_start() {
                exon_ranges.push((exon_start, exon_end));
                maybe_exon_start = None;
            }
        }

        if maybe_exon_start == None {
            maybe_exon_start = Some(rec.get_start());
        }
        maybe_exon_end = Some(rec.get_end());
    }

    if let (Some(exon_start), Some(exon_end)) = (maybe_exon_start, maybe_exon_end) {
        exon_ranges.push((exon_start, exon_end));
    }

    let transcript_start = maybe_transcript_start.unwrap();
    let transcript_end = maybe_transcript_end.unwrap();

    let mut recs_out = Vec::with_capacity(2 + recs.len() * 2);

    let gene_attributes = format!("ID={}", gene_name);
    let gene_rec = GffRecord::new(
        sequence.to_owned(),
        source.to_owned(),
        GffFeature::Gene,
        transcript_start,
        transcript_end,
        None,
        strand,
        None,
        gene_attributes,
    );
    recs_out.push(gene_rec);

    let mrna_attributes = format!("ID={}.1;Parent={}", gene_name, gene_name);
    let mrna_rec = GffRecord::new(
        sequence.to_owned(),
        source.to_owned(),
        GffFeature::MRNA,
        transcript_start,
        transcript_end,
        None,
        strand,
        None,
        mrna_attributes,
    );
    recs_out.push(mrna_rec);

    let mut exon_range_iter = exon_ranges.into_iter();
    let mut current_exon_end = None;
    let mut exon_idx = 1;

    for rec in recs {
        if current_exon_end.is_none() || current_exon_end.unwrap() < rec.get_start() {
            let (exon_start, exon_end) = exon_range_iter.next().unwrap();
            let exon_attributes = format!(
                "ID={}.1.exon.{};Parent={}.1",
                gene_name, exon_idx, gene_name
            );

            let exon_rec = GffRecord::new(
                sequence.to_owned(),
                source.to_owned(),
                GffFeature::Exon,
                exon_start,
                exon_end,
                None,
                strand,
                None,
                exon_attributes,
            );

            recs_out.push(exon_rec);

            exon_idx += 1;
            current_exon_end = Some(exon_end);
        }

        recs_out.push(rec);
    }

    recs_out
}

fn convert_regions_to_gff(
    regions: Vec<HmmStateRegion>,
    sequence: &str,
    source: &str,
    strand: Option<GffStrand>,
    position: usize,
    gene_name: &str,
) -> Vec<GffRecord> {
    // Unlike the numbered Phytozome IDs above, all lines of one CDS, 5' UTR or 3' UTR share one ID,
    // as each is a single discontinuous feature in the GFF3 specification; exons stay numbered.
    let mut coding_offset = 0;

    let mut region_vec = Vec::new();

    for region in regions.iter() {
        let maybe_feature_and_attributes = match region.get_annotation_label() {
            HmmAnnotationLabel::Intergenic => None, //panic!("Intergenic should be removed before now"),
            HmmAnnotationLabel::UTR5 => Some((
                GffFeature::FivePrimeUTR,
                format!("ID={}.1.five_prime_UTR;Parent={}.1", gene_name, gene_name),
            )),
            HmmAnnotationLabel::Start => panic!("Start should be Coding"),
            HmmAnnotationLabel::Coding => Some((
                GffFeature::CDS,
                format!("ID={}.1.CDS;Parent={}.1", gene_name, gene_name),
            )),
            HmmAnnotationLabel::Intron => None,
            HmmAnnotationLabel::Stop => panic!("Stop should be Coding"),
            HmmAnnotationLabel::UTR3 => Some((
                GffFeature::ThreePrimeUTR,
                format!("ID={}.1.three_prime_UTR;Parent={}.1", gene_name, gene_name),
            )),
        };

        if let Some((feature, attributes)) = maybe_feature_and_attributes {
            let start = (region.get_start_pos() + position + 1) as u64;
            let end = (region.get_end_pos() + position) as u64;

            let phase = if region.get_annotation_label() == HmmAnnotationLabel::Coding {
                Some(GffPhase::from(coding_offset))
            } else {
                None
            };

            let rec = GffRecord::new(
                sequence.to_owned(),
                source.to_owned(),
                feature,
                start,
                end,
                None,
                strand,
                phase,
                attributes,
            );

            region_vec.push(rec);

            if region.get_annotation_label() == HmmAnnotationLabel::Coding {
                coding_offset += region.len() as u64;
            }
        }
    }

    region_vec
}

pub fn hmm_solution_to_gff(
    genes: Vec<(Vec<HmmStateRegion>, usize)>,
    species: &str,
    sequence: &str,
    source: &str,
    rev: bool,
    position: usize,
    sequence_length: u64,
    min_coding_length: usize,
    gene_idx: &mut usize,
) -> Vec<GffRecord> {
    //        if genes.len() > 1
    //            { println!("{} genes at {} in {}", genes.len(), position, sequence); }

    let mut all_gff_recs = Vec::new();

    let strand = Some(GffStrand::Forward); // Initially generate everything as forward

    for (gene_regions, coding_length) in genes {
        if coding_length >= min_coding_length {
            let gene_name = format!("{}_{}_{:06}", species, sequence, *gene_idx);
            let gene_gff_recs = convert_regions_to_gff(
                gene_regions,
                sequence,
                source,
                strand,
                position,
                &gene_name,
            );
            let gene_gff_recs =
                generate_gff_aggregate_records(gene_gff_recs, sequence, source, strand, &gene_name);

            all_gff_recs.extend(gene_gff_recs);
            *gene_idx += 1;
        }
    }

    if rev {
        for gff_recs in all_gff_recs.iter_mut() {
            gff_recs.swap_strand(sequence_length);
        }
    }

    all_gff_recs
}
