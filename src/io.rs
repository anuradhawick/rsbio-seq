//! Buffered FASTA/FASTQ parsing inspired by klib's kseq.h:
//! https://github.com/attractivechaos/klib/blob/master/kseq.h
//!
//! Reuse byte buffers, concatenate sequence lines, and consume quality by
//! sequence length. This is an independent Rust implementation with strict
//! format validation rather than kseq's permissive header scanning.
use crate::seq::{SeqFormat, Sequence};
use std::io::{self, BufRead};

pub struct FastxReader<R> {
    reader: R,
    format: SeqFormat,
    line: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    pending_header: bool,
    finished: bool,
}

impl<R: BufRead> FastxReader<R> {
    pub fn new(format: SeqFormat, reader: R) -> Self {
        Self {
            reader,
            format,
            line: Vec::with_capacity(256),
            seq: Vec::with_capacity(256),
            qual: Vec::with_capacity(256),
            pending_header: false,
            finished: false,
        }
    }

    // BufRead::read_until scans buffered chunks for LF and copies slices in
    // bulk. No per-byte I/O, UTF-8 decoding, or per-line allocation is needed.
    fn read_line(&mut self) -> io::Result<bool> {
        self.line.clear();
        if self.reader.read_until(b'\n', &mut self.line)? == 0 {
            return Ok(false);
        }
        if self.line.last() == Some(&b'\n') {
            self.line.pop();
            if self.line.last() == Some(&b'\r') {
                self.line.pop();
            }
        }
        Ok(true)
    }

    fn invalid(message: &str) -> io::Error {
        io::Error::new(io::ErrorKind::InvalidData, message)
    }

    /// Read one record. EOF and errors are terminal; errors never produce a
    /// partially parsed record. Returned strings own their data while scratch
    /// buffers retain their capacity for the next call.
    pub fn read_next(&mut self) -> io::Result<Option<Sequence>> {
        if self.finished {
            return Ok(None);
        }
        let result = self.read_record();
        if !matches!(result, Ok(Some(_))) {
            self.finished = true;
        }
        result
    }

    fn read_record(&mut self) -> io::Result<Option<Sequence>> {
        if self.pending_header {
            self.pending_header = false;
        } else {
            loop {
                if !self.read_line()? {
                    return Ok(None);
                }
                if !self.line.is_empty() {
                    break;
                }
            }
        }
        let marker = match self.format {
            SeqFormat::Fasta => b'>',
            SeqFormat::Fastq => b'@',
        };
        if self.line.first() != Some(&marker) {
            return Err(Self::invalid("Expected a FASTA/FASTQ header"));
        }
        let header = &self.line[1..];
        let split = header
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(header.len());
        if split == 0 {
            return Err(Self::invalid("Sequence identifier is empty"));
        }
        let id = String::from_utf8_lossy(&header[..split]).into_owned();
        let desc = if split < header.len() {
            String::from_utf8_lossy(&header[split + 1..]).into_owned()
        } else {
            String::new()
        };
        self.seq.clear();
        self.qual.clear();
        match self.format {
            SeqFormat::Fasta => {
                while self.read_line()? {
                    if self.line.first() == Some(&b'>') {
                        self.pending_header = true;
                        break;
                    }
                    self.seq.extend_from_slice(&self.line);
                }
            }
            SeqFormat::Fastq => {
                loop {
                    if !self.read_line()? {
                        return Err(Self::invalid("FASTQ record is missing its '+' separator"));
                    }
                    if self.line.first() == Some(&b'+') {
                        break;
                    }
                    if self.line.first() == Some(&b'@') {
                        return Err(Self::invalid("FASTQ record is missing its '+' separator"));
                    }
                    self.seq.extend_from_slice(&self.line);
                }
                // Read at least one quality line, including for empty records.
                loop {
                    if !self.read_line()? {
                        return Err(Self::invalid("FASTQ quality string is truncated"));
                    }
                    self.qual.extend_from_slice(&self.line);
                    if self.qual.len() >= self.seq.len() {
                        break;
                    }
                }
                if self.qual.len() != self.seq.len() {
                    return Err(Self::invalid("FASTQ sequence and quality lengths differ"));
                }
            }
        }
        Ok(Some(Sequence {
            id,
            desc,
            seq: String::from_utf8_lossy(&self.seq).into_owned(),
            qual: String::from_utf8_lossy(&self.qual).into_owned(),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor, Read};

    fn parse(format: SeqFormat, input: &[u8], capacity: usize) -> io::Result<Vec<Sequence>> {
        let mut reader = FastxReader::new(format, BufReader::with_capacity(capacity, input));
        let mut records = Vec::new();
        while let Some(record) = reader.read_next()? {
            records.push(record);
        }
        assert!(reader.read_next()?.is_none());
        Ok(records)
    }

    #[test]
    fn wrapped_records_and_buffer_boundaries() {
        for capacity in [1, 2, 7, 64 * 1024] {
            let fasta = parse(
                SeqFormat::Fasta,
                b"\n>one\tdescription with spaces\r\nAC\r\n\r\nGT\r\n>empty\r\n>last\nN",
                capacity,
            )
            .unwrap();
            assert_eq!(fasta.len(), 3);
            assert_eq!(
                fasta[0],
                Sequence::new(
                    "one".into(),
                    "description with spaces".into(),
                    "ACGT".into(),
                    "".into()
                )
            );
            assert!(fasta[1].seq.is_empty());
            assert_eq!(fasta[2].seq, "N");
            let fastq = parse(
                SeqFormat::Fastq,
                b"@one desc\r\nAC\r\nGT\r\n+one\r\n@+\r\n!I\r\n@two\nN\n+\n#",
                capacity,
            )
            .unwrap();
            assert_eq!(fastq.len(), 2);
            assert_eq!(fastq[0].seq, "ACGT");
            assert_eq!(fastq[0].qual, "@+!I");
            assert_eq!(fastq[1].qual, "#");
        }
    }

    #[test]
    fn long_lines_and_owned_records() {
        let bases = "ACGT".repeat(100_000);
        let input = format!(">long\n{}\n>short\nN\n", bases);
        let records = parse(SeqFormat::Fasta, input.as_bytes(), 17).unwrap();
        assert_eq!(records[0].seq, bases);
        assert_eq!(records[1].seq, "N");
    }

    #[test]
    fn empty_input_and_empty_fastq_record() {
        assert!(parse(SeqFormat::Fasta, b"", 1).unwrap().is_empty());
        let records = parse(SeqFormat::Fastq, b"@empty\n\n+\n\n@next\nA\n+\n!\n", 1).unwrap();
        assert!(records[0].seq.is_empty());
        assert!(records[0].qual.is_empty());
        assert_eq!(records[1].seq, "A");
    }

    #[test]
    fn malformed_records_fail_once() {
        for input in [
            b"bad\n".as_slice(),
            b"@\nA\n+\n!\n",
            b"@x\nAC\n",
            b"@x\nAC\n@next\n",
            b"@x\nAC\n+\n!",
            b"@x\nA\n+\n!!\n",
            b"@x\nA\n+",
        ] {
            let mut reader = FastxReader::new(SeqFormat::Fastq, Cursor::new(input));
            assert_eq!(
                reader.read_next().unwrap_err().kind(),
                io::ErrorKind::InvalidData
            );
            assert!(reader.read_next().unwrap().is_none());
        }
        assert!(parse(SeqFormat::Fasta, b">\nAC\n", 1).is_err());
    }

    #[test]
    fn io_errors_propagate() {
        struct Broken;
        impl Read for Broken {
            fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
                Err(io::Error::other("broken stream"))
            }
        }
        let mut reader = FastxReader::new(SeqFormat::Fasta, BufReader::new(Broken));
        assert_eq!(reader.read_next().unwrap_err().to_string(), "broken stream");
        assert!(reader.read_next().unwrap().is_none());
    }

    #[test]
    fn gzip_members_and_bgzf_are_read_completely() {
        use flate2::{write::GzEncoder, Compression};
        use std::io::Write;
        let mut compressed = Vec::new();
        for record in [b">first\nAC\n".as_slice(), b">second\nGT\n"] {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            encoder.write_all(record).unwrap();
            compressed.extend(encoder.finish().unwrap());
        }
        let input = flate2::read::MultiGzDecoder::new(compressed.as_slice());
        let mut reader = FastxReader::new(SeqFormat::Fasta, BufReader::new(input));
        assert_eq!(reader.read_next().unwrap().unwrap().id, "first");
        assert_eq!(reader.read_next().unwrap().unwrap().id, "second");
        assert!(reader.read_next().unwrap().is_none());
        let records: Vec<_> = crate::reader::Sequences::new(
            SeqFormat::Fasta,
            crate::reader::get_reader("test_data/reads_indexable.fa.gz").unwrap(),
        )
        .collect();
        assert_eq!(records.len(), 7);
        let plain: Vec<_> = crate::reader::Sequences::new(
            SeqFormat::Fasta,
            crate::reader::get_reader("test_data/reads_indexable.fa").unwrap(),
        )
        .collect();
        assert_eq!(records, plain);
    }
}
