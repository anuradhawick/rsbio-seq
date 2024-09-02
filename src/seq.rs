use bio::io::fasta::Records as FastaRecords;
use bio::io::fastq::Records as FastqRecords;
use pyo3::prelude::*;
use std::io::{BufRead, BufReader};

// Record set entries of type R, which implement BufRead trait (stdin/file)
pub enum RecordSet<R: BufRead> {
    Fasta(FastaRecords<BufReader<R>>),
    Fastq(FastqRecords<BufReader<R>>),
}

/// Sequence entry
#[pyclass]
#[derive(FromPyObject)]
pub struct Sequence {
    /// sequence id
    #[pyo3(get, set)]
    pub id: String,
    /// sequence description
    #[pyo3(get, set)]
    pub desc: String,
    /// sequence string
    #[pyo3(get, set)]
    pub seq: String,
    /// sequence quality string (for FASTQ)
    #[pyo3(get, set)]
    pub qual: String,
}

#[derive(Debug, Clone, Copy)]
pub enum SeqFormat {
    Fasta,
    Fastq,
}

impl SeqFormat {
    pub fn get(path: &str) -> Option<SeqFormat> {
        let mut path = path;
        if path.ends_with(".gz") {
            path = path.trim_end_matches(".gz");
        }
        if path.ends_with(".fq") || path.ends_with(".fastq") {
            return Some(SeqFormat::Fastq);
        } else if path.ends_with(".fasta") || path.ends_with(".fa") || path.ends_with(".fna") {
            return Some(SeqFormat::Fasta);
        }
        None
    }
}
