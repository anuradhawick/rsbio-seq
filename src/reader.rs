use crate::io::FastxReader;
use crate::seq::{SeqFormat, Sequence};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

pub struct Sequences<R: BufRead> {
    records: FastxReader<R>,
}

impl<R: BufRead> Sequences<R> {
    pub fn new(format: SeqFormat, reader: R) -> Self {
        Self {
            records: FastxReader::new(format, reader),
        }
    }

    pub fn read_next(&mut self) -> io::Result<Option<Sequence>> {
        self.records.read_next()
    }
}

impl<R: BufRead> Iterator for Sequences<R> {
    type Item = Sequence;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_next().expect("Error reading sequence")
    }
}

pub fn get_reader(path: &str) -> Result<BufReader<Box<dyn Read + Sync + Send>>, String> {
    let file = File::open(path).map_err(|e| format!("Unable to open {}: {}", path, e))?;
    let input: Box<dyn Read + Sync + Send> = if path.ends_with(".gz") {
        // Gzip members include BGZF blocks; decode all members, not just the first.
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(BufReader::with_capacity(64 * 1024, input))
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
        let mut seqs = Sequences::new(SeqFormat::Fastq, reader);
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
        let mut seqs = Sequences::new(SeqFormat::Fastq, reader);
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
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader);
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
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader);
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
        let mut seqs = Sequences::new(SeqFormat::Fasta, reader);
        let record_1 = seqs.next().unwrap();
        assert_eq!("Record_1", record_1.id);
        assert_eq!("ACGTACGTACGT", record_1.seq);
        let finish = seqs.next();
        assert!(finish.is_none());
    }
}
