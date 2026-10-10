use crate::seq::{SeqFormat, Sequence};
use bgzip::{BGZFWriter, Compression as BGZFCompression};
use byteorder::{LittleEndian, WriteBytesExt};
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
    Plain((File, Option<File>)),
    Gzip(GzEncoder<File>),
    BGzip((BGZFWriter<File>, File, File)),
}

pub struct Writer {
    writer: Option<WriterType>,
    format: SeqFormat,
    offset: usize,
}

impl Writer {
    pub fn new(format: SeqFormat, writer: WriterType) -> Self {
        Self {
            writer: Some(writer),
            format,
            offset: 0,
        }
    }

    pub fn close(&mut self) -> Result<(), String> {
        if let Some(writer) = self.writer.take() {
            match writer {
                WriterType::Plain((mut w, fai_file)) => {
                    w.flush().map_err(|e| e.to_string())?;
                    if let Some(mut fai_file) = fai_file {
                        fai_file.flush().map_err(|e| e.to_string())?;
                    };
                    Ok(())
                }
                WriterType::Gzip(mut w) => w.try_finish().map_err(|e| e.to_string()),
                WriterType::BGzip((w, mut fai_file, mut gzi_file)) => {
                    let index = w.close().map_err(|e| e.to_string())?;

                    if let Some(index) = index {
                        gzi_file
                            .write_u64::<LittleEndian>(index.entries().len() as u64)
                            .map_err(|e| e.to_string())?;
                        for entry in index.entries() {
                            gzi_file
                                .write_u64::<LittleEndian>(entry.compressed_offset)
                                .map_err(|e| e.to_string())?;
                            gzi_file
                                .write_u64::<LittleEndian>(entry.uncompressed_offset)
                                .map_err(|e| e.to_string())?;
                        }
                        gzi_file.flush().map_err(|e| e.to_string())?;
                        fai_file.flush().map_err(|e| e.to_string())?;
                    }
                    Ok(())
                }
            }
        } else {
            Err("Writer already closed".to_string())
        }
    }

    #[inline]
    pub fn write(&mut self, seq: Sequence, wrap: Option<u32>) -> Result<(), String> {
        let mut writer = self.writer.as_mut().unwrap();
        let (writer, fai_writer) = match &mut writer {
            WriterType::Gzip(gz) => (gz as &mut dyn Write, None),
            WriterType::Plain((file, fai)) => (
                file as &mut dyn Write,
                fai.as_mut().map(|f| f as &mut dyn Write),
            ),
            WriterType::BGzip((bgz, fai, _)) => {
                (bgz as &mut dyn Write, Some(fai as &mut dyn Write))
            }
        };

        match self.format {
            SeqFormat::Fasta => {
                let seq_str = if let Some(wrap) = wrap {
                    if wrap < 10 {
                        return Err("Wrap value must be at least 10".to_string());
                    }
                    wrap_string_no_whitespace(&seq.seq, wrap as usize)
                } else {
                    seq.seq.clone()
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

                // if index requested write it
                if let Some(fai_writer) = fai_writer {
                    self.offset += seq.id.len() + seq.desc.len() + 3;

                    fai_writer
                        .write(
                            format!(
                                "{}\t{}\t{}\t{}\t{}\n",
                                seq.id,
                                seq.seq.len(),
                                self.offset,
                                wrap.unwrap_or(seq_str.len() as u32),
                                wrap.unwrap_or(seq_str.len() as u32) + 1
                            )
                            .as_bytes(),
                        )
                        .map_err(|e| e.to_string())?;
                    self.offset += seq_str.len() + 1;
                }
            }
            SeqFormat::Fastq => {
                let (seq_str, qual_str) = if let Some(wrap) = wrap {
                    if wrap < 10 {
                        return Err("Wrap value must be at least 10".to_string());
                    }
                    (
                        wrap_string_no_whitespace(&seq.seq, wrap as usize),
                        wrap_string_no_whitespace(&seq.qual, wrap as usize),
                    )
                } else {
                    (seq.seq.clone(), seq.qual.clone())
                };
                let mut buffer = Vec::with_capacity(
                    7 + seq.id.len() + seq.desc.len() + seq_str.len() + qual_str.len(),
                );

                buffer.extend_from_slice(b"@");
                buffer.extend_from_slice(seq.id.as_bytes());
                buffer.extend_from_slice(b" ");
                buffer.extend_from_slice(seq.desc.as_bytes());
                buffer.extend_from_slice(b"\n");
                buffer.extend_from_slice(seq_str.as_bytes());
                buffer.extend_from_slice(b"\n+\n");
                buffer.extend_from_slice(qual_str.as_bytes());
                buffer.extend_from_slice(b"\n");
                writer.write_all(&buffer).map_err(|e| e.to_string())?;

                // if index requested write it
                if let Some(fai_writer) = fai_writer {
                    self.offset += seq.id.len() + seq.desc.len() + 3; // 3 = > "@", " ", and "\n"
                    let qual_offset = self.offset + seq_str.len() + 3; // 3 = > "\n+\n"

                    fai_writer
                        .write(
                            format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\n",
                                seq.id,
                                seq.seq.len(),
                                self.offset,
                                wrap.unwrap_or(seq_str.len() as u32),
                                wrap.unwrap_or(seq_str.len() as u32) + 1,
                                qual_offset,
                            )
                            .as_bytes(),
                        )
                        .map_err(|e| e.to_string())?;
                    self.offset += seq_str.len() + 3 + qual_str.len() + 1;
                }
            }
        }
        Ok(())
    }
}

pub fn get_writer(path: &str, index: bool) -> Result<WriterType, String> {
    let is_zip = path.ends_with(".gz");
    let file = File::create(path).map_err(|_| format!("Unable to open: {}", path))?;
    if is_zip {
        if index {
            let encoder = BGZFWriter::new(file, BGZFCompression::default());
            let fai_file = File::create(format!("{}.fai", path))
                .map_err(|_| format!("Unable to open: {}", path))?;
            let gzi_file = File::create(format!("{}.gzi", path))
                .map_err(|_| format!("Unable to open: {}", path))?;

            Ok(WriterType::BGzip((encoder, fai_file, gzi_file)))
        } else {
            let encoder = GzEncoder::new(file, Compression::best());
            Ok(WriterType::Gzip(encoder))
        }
    } else if index {
        let fai_file = File::create(format!("{}.fai", path))
            .map_err(|_| format!("Unable to open: {}", path))?;
        Ok(WriterType::Plain((file, Some(fai_file))))
    } else {
        Ok(WriterType::Plain((file, None)))
    }
}

#[cfg(test)]
mod writer_tests {
    use super::{get_writer, Writer};
    use crate::seq::{SeqFormat, Sequence};
    use flate2::read::{GzDecoder, MultiGzDecoder};
    use std::{fs::File, io::Read};

    const PATH_FQ: &str = "test_data/out.fq";
    const PATH_FA: &str = "test_data/out.fa";
    const PATH_FA_IX: &str = "test_data/out_ix.fa";
    const PATH_FQ_GZ: &str = "test_data/out.fq.gz";
    const PATH_FA_GZ: &str = "test_data/out.fa.gz";
    const PATH_FA_GZ_IX: &str = "test_data/out_ix.fa.gz";

    fn assert_indexed_fastq_output(path: &str, wrap: Option<u32>) {
        let mut writer = Writer::new(SeqFormat::Fastq, get_writer(path, true).unwrap());
        writer
            .write(
                Sequence {
                    id: "rec_1".into(),
                    desc: "desc 1".into(),
                    seq: "ACGTACGTACGT".into(),
                    qual: "!5I@+!5I@+!5".into(),
                },
                wrap,
            )
            .unwrap();
        writer
            .write(
                Sequence {
                    id: "rec_2".into(),
                    desc: "desc 2".into(),
                    seq: "TGCATGCATGCATGCATGCA".into(),
                    qual: "I".repeat(20),
                },
                wrap,
            )
            .unwrap();
        writer.close().unwrap();

        let compressed = path.ends_with(".gz");
        let mut output = String::new();
        let file = File::open(path).unwrap();
        if compressed {
            MultiGzDecoder::new(file)
                .read_to_string(&mut output)
                .unwrap();
            let header = std::fs::read(path).unwrap();
            assert_eq!(&header[..3], &[0x1f, 0x8b, 0x08]);
            assert_eq!(&header[12..14], b"BC");
            // These small records fit in one BGZF block; no extra offsets.
            assert_eq!(std::fs::read(format!("{path}.gzi")).unwrap(), [0; 8]);
        } else {
            let mut file = file;
            file.read_to_string(&mut output).unwrap();
        }

        let (expected_output, expected_index) = if wrap.is_some() {
            (
                "@rec_1 desc 1\nACGTACGTAC\nGT\n+\n!5I@+!5I@+\n!5\n\
@rec_2 desc 2\nTGCATGCATG\nCATGCATGCA\n+\nIIIIIIIIII\nIIIIIIIIII\n",
                "rec_1\t12\t14\t10\t11\t30\nrec_2\t20\t58\t10\t11\t82\n",
            )
        } else {
            (
                "@rec_1 desc 1\nACGTACGTACGT\n+\n!5I@+!5I@+!5\n\
@rec_2 desc 2\nTGCATGCATGCATGCATGCA\n+\nIIIIIIIIIIIIIIIIIIII\n",
                "rec_1\t12\t14\t12\t13\t29\nrec_2\t20\t56\t20\t21\t79\n",
            )
        };
        assert_eq!(output, expected_output);
        assert_eq!(
            std::fs::read_to_string(format!("{path}.fai")).unwrap(),
            expected_index
        );
    }

    #[test]
    fn write_fq_file_indexed_test() {
        assert_indexed_fastq_output("test_data/out_ix.fq", Some(10));
    }

    #[test]
    fn write_fq_gz_file_indexed_test() {
        assert_indexed_fastq_output("test_data/out_ix.fq.gz", Some(10));
    }

    #[test]
    fn write_fq_file_indexed_unwrapped_test() {
        assert_indexed_fastq_output("test_data/out_ix_unwrapped.fq", None);
    }

    #[test]
    fn write_fq_gz_file_indexed_unwrapped_test() {
        assert_indexed_fastq_output("test_data/out_ix_unwrapped.fq.gz", None);
    }

    #[test]
    fn write_fa_file_test() {
        // write
        let writer = get_writer(PATH_FA, false).unwrap();
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
    fn write_fa_file_indexed_test() {
        // write
        let writer = get_writer(PATH_FA_IX, true).unwrap();
        let format = SeqFormat::get(PATH_FA_IX).unwrap();
        let seq_1 = Sequence {
            id: "rec_1".into(),
            desc: "desc 1".into(),
            seq: "ACGTCCGTACGTCCGTTTCGACGTCCGTACGTCCGTTTCGACGTCCGTACGTCCGTTTCG".into(),
            qual: "".into(),
        };
        let seq_2 = Sequence {
            id: "rec_2".into(),
            desc: "desc 2".into(),
            seq: "CCGTACGTCCGTCCGTACGTCCGTGTCCGTCCGTACGTCCGTGTCCGTCCGTACGTCCGTGTCCGTCCGTACGTCCGT"
                .into(),
            qual: "".into(),
        };
        let mut seq_writer = Writer::new(format, writer);
        seq_writer.write(seq_1, Some(10)).unwrap();
        seq_writer.write(seq_2, Some(10)).unwrap();
        seq_writer.close().unwrap();
        // read and validate
        let mut buf = Vec::new();
        File::open(format!("{}.fai", PATH_FA_IX))
            .unwrap()
            .read_to_end(&mut buf)
            .unwrap();
        assert_eq!(
            String::from_utf8_lossy(&buf).to_string(),
            "rec_1\t60\t14\t10\t11\nrec_2\t78\t94\t10\t11\n"
        );
    }

    #[test]
    fn write_fa_gz_file_indexed_test() {
        // write
        let writer = get_writer(PATH_FA_GZ_IX, true).unwrap();
        let format = SeqFormat::get(PATH_FA_GZ_IX).unwrap();
        let seq_1 = Sequence {
            id: "rec_1".into(),
            desc: "desc 1".into(),
            seq: "ACGTCCGTACGTCCGTTTCGACGTCCGTACGTCCGTTTCGACGTCCGTACGTCCGTTTCG".into(),
            qual: "".into(),
        };
        let seq_2 = Sequence {
            id: "rec_2".into(),
            desc: "desc 2".into(),
            seq: "CCGTACGTCCGTCCGTACGTCCGTGTCCGTCCGTACGTCCGTGTCCGTCCGTACGTCCGTGTCCGTCCGTACGTCCGT"
                .into(),
            qual: "".into(),
        };
        let mut seq_writer = Writer::new(format, writer);
        seq_writer.write(seq_1, Some(10)).unwrap();
        seq_writer.write(seq_2, Some(10)).unwrap();
        seq_writer.close().unwrap();
        // read and validate
        let mut buf = Vec::new();
        File::open(format!("{}.fai", PATH_FA_GZ_IX))
            .unwrap()
            .read_to_end(&mut buf)
            .unwrap();
        assert_eq!(
            String::from_utf8_lossy(&buf).to_string(),
            "rec_1\t60\t14\t10\t11\nrec_2\t78\t94\t10\t11\n"
        );
    }

    #[test]
    fn write_fq_file_test() {
        // write
        let writer = get_writer(PATH_FQ, false).unwrap();
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
        let writer = get_writer(PATH_FA_GZ, false).unwrap();
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
        let writer = get_writer(PATH_FQ_GZ, false).unwrap();
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
