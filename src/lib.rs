mod seq;
use pyo3::prelude::*;
use seq::{get_reader, SeqFormat, Sequence, Sequences};
use std::io::Read;

/// Sequence reader
#[pyclass]
pub struct SeqIO {
    records: Sequences<std::io::BufReader<Box<dyn Read + Send + Sync>>>,
}

#[pymethods]
impl SeqIO {
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

/// Sequence reader for rust
#[pymodule]
fn rsbio_seq(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Sequence>()?;
    m.add_class::<SeqIO>()?;
    Ok(())
}
