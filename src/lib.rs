mod reader;
mod seq;
mod writer;
use pyo3::{exceptions::PyIOError, prelude::*};
use reader::{get_reader, Sequences};
use seq::{SeqFormat, Sequence};
use std::io::Read;
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
    pub fn new(path: String) -> PyResult<Self> {
        let reader = get_reader(&path).map_err(PyIOError::new_err)?;
        let format = SeqFormat::get(&path).map_err(PyIOError::new_err)?;

        Ok(Self {
            records: Sequences::new(format, reader),
        })
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
    writer: Writer,
}

#[pymethods]
impl SeqWriter {
    /// Initialise a sequence writer for a file in some destination
    #[new]
    #[pyo3(signature = (path))]
    pub fn new(path: String) -> PyResult<Self> {
        let writer = get_writer(&path).map_err(PyIOError::new_err)?;
        let format = SeqFormat::get(&path).map_err(PyIOError::new_err)?;

        Ok(Self {
            writer: Writer::new(format, writer),
        })
    }

    #[pyo3(signature = (seq, wrap=None))]
    pub fn write(&mut self, seq: Sequence, wrap: Option<u32>) -> PyResult<()> {
        self.writer.write(seq, wrap).map_err(PyIOError::new_err)
    }

    #[pyo3(signature = ())]
    pub fn close(&mut self) -> PyResult<()> {
        self.writer.close().map_err(PyIOError::new_err)
    }
}

#[pyfunction]
pub fn phred_to_ascii(scores: Vec<u8>) -> PyResult<String> {
    Ok(scores
        .iter()
        .map(|&score| (score + 33) as char) // Convert Phred score to ASCII
        .collect())
}

#[pyfunction]
pub fn ascii_to_phred(qual: String) -> PyResult<Vec<u8>> {
    Ok(qual
        .chars()
        .map(|c| (c as u8).saturating_sub(33)) // Convert ASCII to Phred score
        .collect())
}

/// Sequence reader for rust
#[pymodule]
fn rsbio_seq(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Sequence>()?;
    m.add_class::<SeqReader>()?;
    m.add_class::<SeqWriter>()?;
    m.add_function(wrap_pyfunction!(phred_to_ascii, m)?)?;
    m.add_function(wrap_pyfunction!(ascii_to_phred, m)?)?;
    Ok(())
}
