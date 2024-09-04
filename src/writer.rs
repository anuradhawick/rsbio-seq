use crate::seq::{SeqFormat, Sequence};
use flate2::{write::GzEncoder, Compression};
use std::{fs::File, io::Write};

#[inline]
fn wrap_string_no_whitespace(s: &str, width: usize) -> String {
    let mut result = String::with_capacity(s.len() + s.len() / width);
    let mut i = 0;

    while i < s.len() {
        let end = if i + width < s.len() {
            i + width
        } else {
            s.len()
        };
        result.push_str(&s[i..end]);
        result.push('\n');
        i += width;
    }

    // Remove the last newline
    if result.ends_with('\n') {
        result.pop();
    }

    result
}

pub enum WriterType {
    Plain(File),
    Gzip(GzEncoder<File>),
}

pub struct Writer {
    writer: WriterType,
    format: SeqFormat,
}

impl Writer {
    pub fn new(format: SeqFormat, writer: WriterType) -> Self {
        Self { writer, format }
    }

    pub fn close(&mut self) -> Result<(), String> {
        match &mut self.writer {
            WriterType::Plain(w) => w.flush().map_err(|e| e.to_string()),
            WriterType::Gzip(w) => w.try_finish().map_err(|e| e.to_string()),
        }
    }

    #[inline]
    pub fn write(&mut self, seq: Sequence, wrap: Option<u32>) -> Result<(), String> {
        let writer = match &mut self.writer {
            WriterType::Gzip(gz) => gz as &mut dyn Write,
            WriterType::Plain(file) => file as &mut dyn Write,
        };

        match self.format {
            SeqFormat::Fasta => {
                let seq_str = if let Some(wrap) = wrap {
                    if wrap < 10 {
                        return Err("Wrap value must be at least 10".to_string());
                    }
                    wrap_string_no_whitespace(&seq.seq, wrap as usize)
                } else {
                    seq.seq
                };
                let mut buffer =
                    Vec::with_capacity(4 + seq.id.len() + seq.desc.len() + seq_str.len());

                buffer.extend_from_slice(b">");
                buffer.extend_from_slice(seq.id.as_bytes());
                buffer.extend_from_slice(b" ");
                buffer.extend_from_slice(seq.desc.as_bytes());
                buffer.extend_from_slice(b"\n");
                buffer.extend_from_slice(seq_str.as_bytes());
                buffer.extend_from_slice(b"\n");
                writer.write_all(&buffer).map_err(|e| e.to_string())?;
            }
            SeqFormat::Fastq => {
                let mut buffer = Vec::with_capacity(
                    7 + seq.id.len() + seq.desc.len() + seq.seq.len() + seq.qual.len(),
                );

                buffer.extend_from_slice(b"@");
                buffer.extend_from_slice(seq.id.as_bytes());
                buffer.extend_from_slice(b" ");
                buffer.extend_from_slice(seq.desc.as_bytes());
                buffer.extend_from_slice(b"\n");
                buffer.extend_from_slice(seq.seq.as_bytes());
                buffer.extend_from_slice(b"\n+\n");
                buffer.extend_from_slice(seq.qual.as_bytes());
                buffer.extend_from_slice(b"\n");
                writer.write_all(&buffer).map_err(|e| e.to_string())?;
            }
        }
        Ok(())
    }
}

pub fn get_writer(path: &str) -> Result<WriterType, String> {
    let is_zip = path.ends_with(".gz");
    let file = File::create(path).map_err(|_| format!("Unable to open: {}", path))?;
    if is_zip {
        let encoder = GzEncoder::new(file, Compression::best());
        Ok(WriterType::Gzip(encoder))
    } else {
        Ok(WriterType::Plain(file))
    }
}

#[cfg(test)]
mod writer_tests {
    use super::{get_writer, Writer};
    use crate::seq::{SeqFormat, Sequence};
    use flate2::read::GzDecoder;
    use std::{fs::File, io::Read};

    const PATH_FQ: &str = "test_data/out.fq";
    const PATH_FA: &str = "test_data/out.fa";
    const PATH_FQ_GZ: &str = "test_data/out.fq.gz";
    const PATH_FA_GZ: &str = "test_data/out.fa.gz";

    #[test]
    fn write_fa_file_test() {
        // write
        let writer = get_writer(PATH_FA).unwrap();
        let format = SeqFormat::get(PATH_FA).unwrap();
        let seq = Sequence {
            id: "rec_1".into(),
            desc: "desc 1".into(),
            seq: "ACGTCCGT".into(),
            qual: "".into(),
        };
        let mut seq_writer = Writer::new(format, writer);
        seq_writer.write(seq, None).unwrap();
        seq_writer.close().unwrap();
        // read and validate
        let mut buf = Vec::new();
        File::open(PATH_FA).unwrap().read_to_end(&mut buf).unwrap();
        assert_eq!(
            String::from_utf8_lossy(&buf).to_string(),
            ">rec_1 desc 1\nACGTCCGT\n"
        );
    }

    #[test]
    fn write_fq_file_test() {
        // write
        let writer = get_writer(PATH_FQ).unwrap();
        let format = SeqFormat::get(PATH_FQ).unwrap();
        let seq = Sequence {
            id: "rec_1".into(),
            desc: "desc 1".into(),
            seq: "ACGTCCGT".into(),
            qual: "IIIIIIII".into(),
        };
        let mut seq_writer = Writer::new(format, writer);
        seq_writer.write(seq, None).unwrap();
        seq_writer.close().unwrap();
        // read and validate
        let mut buf = Vec::new();
        File::open(PATH_FQ).unwrap().read_to_end(&mut buf).unwrap();
        assert_eq!(
            String::from_utf8_lossy(&buf).to_string(),
            "@rec_1 desc 1\nACGTCCGT\n+\nIIIIIIII\n"
        );
    }

    #[test]
    fn write_fa_gz_file_test() {
        // write
        let writer = get_writer(PATH_FA_GZ).unwrap();
        let format = SeqFormat::get(PATH_FA_GZ).unwrap();
        let seq = Sequence {
            id: "rec_1".into(),
            desc: "desc 1".into(),
            seq: "ACGTCCGT".into(),
            qual: "".into(),
        };
        let mut seq_writer = Writer::new(format, writer);
        seq_writer.write(seq, None).unwrap();
        seq_writer.close().unwrap();
        // read and validate
        let mut buf = Vec::new();
        File::open(PATH_FA_GZ)
            .unwrap()
            .read_to_end(&mut buf)
            .unwrap();
        let mut ubuf = String::new();
        let mut gz = GzDecoder::new(&buf[..]);
        gz.read_to_string(&mut ubuf).unwrap();
        assert_eq!(ubuf, ">rec_1 desc 1\nACGTCCGT\n");
    }

    #[test]
    fn write_fq_gz_file_test() {
        // write
        let writer = get_writer(PATH_FQ_GZ).unwrap();
        let format = SeqFormat::get(PATH_FQ_GZ).unwrap();
        let seq = Sequence {
            id: "rec_1".into(),
            desc: "desc 1".into(),
            seq: "ACGTCCGT".into(),
            qual: "IIIIIIII".into(),
        };
        let mut seq_writer = Writer::new(format, writer);
        seq_writer.write(seq, None).unwrap();
        seq_writer.close().unwrap();
        // read and validate
        let mut buf = Vec::new();
        File::open(PATH_FQ_GZ)
            .unwrap()
            .read_to_end(&mut buf)
            .unwrap();
        let mut ubuf = String::new();
        let mut gz = GzDecoder::new(&buf[..]);
        gz.read_to_string(&mut ubuf).unwrap();
        assert_eq!(ubuf, "@rec_1 desc 1\nACGTCCGT\n+\nIIIIIIII\n");
    }
}
