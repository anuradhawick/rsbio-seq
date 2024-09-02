use bio::io::fasta::Reader as FastaReader;
use bio::io::fastq::Reader as FastqReader;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use crate::seq::{RecordSet, SeqFormat, Sequence};

pub struct Sequences<R: BufRead> {
    pub records: RecordSet<R>,
}

impl<R: BufRead> Sequences<R> {
    pub fn new(format: SeqFormat, reader: R) -> Result<Self, String> {
        match format {
            SeqFormat::Fastq => {
                let fastq_reader = FastqReader::new(reader);
                Ok(Sequences {
                    records: RecordSet::Fastq(fastq_reader.records()),
                })
            }
            SeqFormat::Fasta => {
                let fasta_reader = FastaReader::new(reader);
                Ok(Sequences {
                    records: RecordSet::Fasta(fasta_reader.records()),
                })
            }
        }
    }
}

impl<R: BufRead> Iterator for Sequences<R> {
    type Item = Sequence;

    fn next(&mut self) -> Option<Self::Item> {
        // records do not have a common trait to get id and seq, we can create one
        // but this looks simpler for the time being
        match self.records {
            RecordSet::Fastq(ref mut records) => {
                let next_record = records.next();
                if let Some(record) = next_record {
                    let record = record.unwrap();
                    return Some(Sequence {
                        id: record.id().to_string(),
                        desc: record.desc().unwrap_or("").to_string(),
                        seq: String::from_utf8_lossy(record.seq()).to_string(),
                        qual: String::from_utf8_lossy(record.qual()).to_string(),
                    });
                }
                None
            }
            RecordSet::Fasta(ref mut records) => {
                let next_record = records.next();
                if let Some(record) = next_record {
                    let record = record.unwrap();
                    return Some(Sequence {
                        id: record.id().to_string(),
                        desc: record.desc().unwrap_or("").to_string(),
                        seq: String::from_utf8_lossy(record.seq()).to_string(),
                        qual: "".into(),
                    });
                }
                None
            }
        }
    }

    fn count(self) -> usize
    where
        Self: Sized,
    {
        unimplemented!("Count cannot be performed without always having a rewindable input stream, stdin is not!");
    }
}

pub fn get_reader(path: &str) -> Result<BufReader<Box<dyn Read + Sync + Send>>, String> {
    let is_zip = path.ends_with(".gz");
    let file = File::open(path).map_err(|_| format!("Unable to open: {}", path))?;
    if is_zip {
        let decoder = flate2::read::GzDecoder::new(file);
        Ok(BufReader::new(Box::new(decoder)))
    } else {
        Ok(BufReader::new(Box::new(file)))
    }
}

#[cfg(test)]
mod reader_tests {
    use super::*;
    const PATH_FQ: &str = "test_data/reads.fq";
    const PATH_FA: &str = "test_data/reads.fa";
    const PATH_FQ_GZ: &str = "test_data/reads.fq.gz";
    const PATH_FA_GZ: &str = "test_data/reads.fa.gz";

    #[test]
    fn load_fq_file_test() {
        let reader = get_reader(PATH_FQ).unwrap();
        let mut seqs = Sequences::new(SeqFormat::Fastq, reader).unwrap();
        let record_1 = seqs.next().unwrap();
        assert_eq!("Read_1", record_1.id);
        assert_eq!(
            "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
            record_1.seq
        );
        let record_2 = seqs.next().unwrap();
        assert_eq!("Read_2", record_2.id);
        assert_eq!(
            "GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT",
            record_2.seq
        );
        let finish = seqs.next();
        assert!(finish.is_none());
    }

    #[test]
    fn load_fq_gz_file_test() {
        let reader = get_reader(PATH_FQ_GZ).unwrap();
        let mut seqs = Sequences::new(SeqFormat::Fastq, reader).unwrap();
        let record_1 = seqs.next().unwrap();
        assert_eq!("Read_1", record_1.id);
        assert_eq!(
            "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
            record_1.seq
        );
        let record_2 = seqs.next().unwrap();
        assert_eq!("Read_2", record_2.id);
        assert_eq!(
            "GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT",
            record_2.seq
        );
        let finish = seqs.next();
        assert!(finish.is_none());
    }

    #[test]
    fn load_fa_file_test() {
        let reader = get_reader(PATH_FA).unwrap();
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader).unwrap();
        let record_1 = seqs.next().unwrap();
        assert_eq!("Record_1", record_1.id);
        assert_eq!(
            "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
            record_1.seq
        );
        let record_2 = seqs.next().unwrap();
        assert_eq!("Record_2", record_2.id);
        assert_eq!(
            "GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT",
            record_2.seq
        );
        let finish = seqs.next();
        assert!(finish.is_none());
    }

    #[test]
    fn load_fa_gz_file_test() {
        let reader = get_reader(PATH_FA_GZ).unwrap();
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader).unwrap();
        let record_1 = seqs.next().unwrap();
        assert_eq!("Record_1", record_1.id);
        assert_eq!(
            "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA",
            record_1.seq
        );
        let record_2 = seqs.next().unwrap();
        assert_eq!("Record_2", record_2.id);
        assert_eq!(
            "GTTCAGGGATACGACGTTTGTATTTTAAGAATCTGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCAT",
            record_2.seq
        );
        let finish = seqs.next();
        assert!(finish.is_none());
    }

    #[test]
    fn parser_test() {
        let input = ">Record_1\nACGTACGTACGT";
        let reader = BufReader::new(input.as_bytes());
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader).unwrap();
        let record_1 = seqs.next().unwrap();
        assert_eq!("Record_1", record_1.id);
        assert_eq!("ACGTACGTACGT", record_1.seq);
        let finish = seqs.next();
        assert!(finish.is_none());
    }
}
