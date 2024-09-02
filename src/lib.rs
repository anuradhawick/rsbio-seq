mod reader;
mod seq;
mod writer;
use pyo3::prelude::*;
use reader::{get_reader, Sequences};
use seq::{SeqFormat, Sequence};
use std::io::{Read, Write};
use writer::{get_writer, Writer};

/// Sequence reader
#[pyclass]
pub struct SeqReader {
    records: Sequences<std::io::BufReader<Box<dyn Read + Send + Sync>>>,
}

#[pymethods]
impl SeqReader {
    /// Initialise a sequence reader for a file in some destination
    #[new]
    #[pyo3(signature = (path))]
    pub fn new(path: String) -> Self {
        let reader = get_reader(&path).unwrap();
        let format = SeqFormat::get(&path).expect("Unable to detect file format");

        Self {
            records: Sequences::new(format, reader).unwrap(),
        }
    }

    /// Iterator object
    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Iterate through sequences
    pub fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Sequence> {
        slf.records.next()
    }
}

/// Sequence reader
#[pyclass]
pub struct SeqWriter {
    // records: Sequences<std::io::BufReader<Box<dyn Read + Send + Sync>>>,
    writer: Writer<Box<dyn Write + Send + Sync>>,
}

#[pymethods]
impl SeqWriter {
    /// Initialise a sequence reader for a file in some destination
    #[new]
    #[pyo3(signature = (path))]
    pub fn new(path: String) -> Self {
        let writer = get_writer(&path).unwrap();
        let format = SeqFormat::get(&path).expect("Unable to detect file format");

        Self {
            writer: Writer::new(format, Box::new(writer)),
        }
    }

    fn test(&self) {}

    #[pyo3(signature = (seq, wrap=None))]
    pub fn write(&mut self, seq: Sequence, wrap: Option<u32>) {
        self.writer.write(seq, wrap).unwrap();
    }
}

/// Sequence reader for rust
#[pymodule]
fn rsbio_seq(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Sequence>()?;
    m.add_class::<SeqReader>()?;
    Ok(())
}
