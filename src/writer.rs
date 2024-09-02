use crate::seq::{SeqFormat, Sequence};
use flate2::Compression;
use std::{
    fs::File,
    io::{BufWriter, Error, Write},
};
pub struct Writer<W: Write> {
    writer: W,
    format: SeqFormat,
}

impl<W: Write> Writer<W> {
    pub fn new(format: SeqFormat, writer: W) -> Self {
        Self { writer, format }
    }

    pub fn write(&mut self, seq: Sequence, wrap: Option<u32>) -> Result<(), Error> {
        match self.format {
            SeqFormat::Fasta => {
                self.writer
                    .write_fmt(format_args!(">{} {}\n{}\n", seq.id, seq.desc, seq.seq))?;
            }
            SeqFormat::Fastq => {
                self.writer.write_fmt(format_args!(
                    "@{} {}\n{}\n+\n{}\n",
                    seq.id, seq.desc, seq.seq, seq.qual
                ))?;
            }
        }
        Ok(())
    }
}

pub fn get_writer(path: &str) -> Result<BufWriter<Box<dyn Write + Sync + Send>>, String> {
    let is_zip = path.ends_with(".gz");
    let file = File::create(path).map_err(|_| format!("Unable to open: {}", path))?;
    if is_zip {
        println!("GZIP");
        let encoder = flate2::write::GzEncoder::new(file, Compression::best());
        Ok(BufWriter::new(Box::new(encoder)))
    } else {
        Ok(BufWriter::new(Box::new(file)))
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
        drop(seq_writer);
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
        drop(seq_writer);
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
        drop(seq_writer);
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
        drop(seq_writer);
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
