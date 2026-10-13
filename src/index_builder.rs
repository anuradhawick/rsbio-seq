//! Build samtools-compatible FAI and GZI indexes without rewriting the input.
use crate::seq::SeqFormat;
use byteorder::{LittleEndian, WriteBytesExt};
use flate2::read::GzDecoder;
use std::{
    collections::HashSet,
    fs::File,
    io::{self, BufRead, BufReader, Read, Write},
};

fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

/// Decode and validate one BGZF member at a time, recording physical block starts.
struct BgzfReader {
    file: File,
    block: Vec<u8>,
    position: usize,
    compressed: u64,
    uncompressed: u64,
    entries: Vec<(u64, u64)>,
}

impl Read for BgzfReader {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        while self.position == self.block.len() {
            let mut header = [0; 12];
            if self.file.read(&mut header[..1])? == 0 {
                return Ok(0);
            }
            self.file.read_exact(&mut header[1..])?;
            if header[..4] != [31, 139, 8, 4] {
                return Err(invalid(
                    "Indexed compressed input must be BGZF, not ordinary gzip",
                ));
            }
            let extra_len = u16::from_le_bytes([header[10], header[11]]) as usize;
            let mut member = header.to_vec();
            member.resize(12 + extra_len, 0);
            self.file.read_exact(&mut member[12..])?;
            let mut size = None;
            let mut p = 12;
            while p + 4 <= member.len() {
                let n = u16::from_le_bytes([member[p + 2], member[p + 3]]) as usize;
                if p + 4 + n > member.len() {
                    return Err(invalid("Invalid BGZF extra field"));
                }
                if &member[p..p + 2] == b"BC" && n == 2 {
                    size = Some(u16::from_le_bytes([member[p + 4], member[p + 5]]) as usize + 1);
                }
                p += 4 + n;
            }
            let size = size.ok_or_else(|| invalid("Gzip input is not BGZF (missing BC field)"))?;
            if p != member.len() || size < member.len() + 8 {
                return Err(invalid("Invalid BGZF block size"));
            }
            let header_len = member.len();
            member.resize(size, 0);
            self.file.read_exact(&mut member[header_len..])?;
            self.block.clear();
            GzDecoder::new(member.as_slice())
                .take(65537)
                .read_to_end(&mut self.block)?;
            if self.block.len() > 65536 {
                return Err(invalid("BGZF block exceeds 64 KiB"));
            }
            if self.compressed != 0 && self.uncompressed != 0 && !self.block.is_empty() {
                self.entries.push((self.compressed, self.uncompressed));
            }
            self.compressed += size as u64;
            self.uncompressed += self.block.len() as u64;
            self.position = 0;
        }
        let n = out.len().min(self.block.len() - self.position);
        out[..n].copy_from_slice(&self.block[self.position..self.position + n]);
        self.position += n;
        Ok(n)
    }
}

struct Lines<R> {
    reader: R,
    bytes: Vec<u8>,
    offset: u64,
    width: usize,
}
impl<R: BufRead> Lines<R> {
    fn next(&mut self) -> io::Result<bool> {
        self.bytes.clear();
        self.width = self.reader.read_until(b'\n', &mut self.bytes)?;
        self.offset += self.width as u64;
        if self.bytes.last() == Some(&b'\n') {
            self.bytes.pop();
            if self.bytes.last() == Some(&b'\r') {
                self.bytes.pop();
            }
        }
        Ok(self.width != 0)
    }
}

#[derive(Default)]
struct Geometry {
    length: u64,
    bases: usize,
    width: usize,
    previous: Option<(usize, usize)>,
}
impl Geometry {
    fn add(&mut self, bases: usize, width: usize) -> io::Result<()> {
        if bases == 0 {
            return Err(invalid(
                "FAI indexing requires nonempty sequence/quality lines",
            ));
        }
        if let Some(previous) = self.previous {
            if previous != (self.bases, self.width) || bases > self.bases {
                return Err(invalid(
                    "FAI indexing requires uniform wrapping except on the final line",
                ));
            }
            let newline = width - bases;
            if newline != 0 && newline != self.width - self.bases {
                return Err(invalid(
                    "FAI indexing requires consistent line endings within a record",
                ));
            }
        } else {
            self.bases = bases;
            self.width = width;
        }
        self.previous = Some((bases, width));
        self.length += bases as u64;
        Ok(())
    }
}

fn fai<R: BufRead>(format: SeqFormat, reader: R) -> io::Result<(Vec<u8>, R)> {
    let mut lines = Lines {
        reader,
        bytes: Vec::new(),
        offset: 0,
        width: 0,
    };
    let mut output = Vec::new();
    let mut names = HashSet::new();
    let mut available = lines.next()?;
    while available {
        let marker = match format {
            SeqFormat::Fasta => b'>',
            SeqFormat::Fastq => b'@',
        };
        if lines.bytes.first() != Some(&marker) {
            return Err(invalid("Expected a FASTA/FASTQ header"));
        }
        let name = lines.bytes[1..]
            .split(u8::is_ascii_whitespace)
            .next()
            .unwrap_or_default();
        if name.is_empty() {
            return Err(invalid("Sequence identifier is empty"));
        }
        let name = std::str::from_utf8(name)
            .map_err(|_| invalid("Sequence identifier is not UTF-8"))?
            .to_owned();
        if !names.insert(name.clone()) {
            return Err(invalid("Duplicate sequence identifier"));
        }
        let offset = lines.offset;
        let mut seq = Geometry::default();
        available = lines.next()?;
        while available {
            if matches!(format, SeqFormat::Fasta) && lines.bytes.first() == Some(&b'>') {
                break;
            }
            if matches!(format, SeqFormat::Fastq) && lines.bytes.first() == Some(&b'+') {
                break;
            }
            if lines.bytes.iter().any(u8::is_ascii_whitespace) || lines.bytes.first() == Some(&b'@')
            {
                return Err(invalid("Invalid sequence line"));
            }
            seq.add(lines.bytes.len(), lines.width)?;
            available = lines.next()?;
        }
        if seq.length == 0 {
            return Err(invalid("FAI indexing requires nonempty sequences"));
        }
        write!(
            output,
            "{}\t{}\t{}\t{}\t{}",
            name, seq.length, offset, seq.bases, seq.width
        )?;
        if matches!(format, SeqFormat::Fastq) {
            if !available || lines.bytes.first() != Some(&b'+') {
                return Err(invalid("FASTQ record is missing its '+' separator"));
            }
            let qual_offset = lines.offset;
            let mut qual = Geometry::default();
            while qual.length < seq.length {
                if !lines.next()? {
                    return Err(invalid("FASTQ quality string is truncated"));
                }
                qual.add(lines.bytes.len(), lines.width)?;
            }
            if qual.length != seq.length
                || qual.bases != seq.bases
                || (qual.length > qual.bases as u64 && qual.width != seq.width)
                || (qual.width > qual.bases && qual.width - qual.bases != seq.width - seq.bases)
            {
                return Err(invalid(
                    "FASTQ quality must match sequence length and wrapping",
                ));
            }
            write!(output, "\t{}", qual_offset)?;
            available = lines.next()?;
        }
        writeln!(output)?;
    }
    Ok((output, lines.reader))
}

/// Create `<path>.fai` and, for BGZF `.gz` input, `<path>.gzi`.
/// Validate all input before opening output files. Ordinary gzip is unsupported.
pub fn build_index(path: &str) -> io::Result<()> {
    let format = SeqFormat::get(path).map_err(|e| invalid(&e))?;
    let file = File::open(path)?;
    let (index, entries) = if path.ends_with(".gz") {
        let reader = BgzfReader {
            file,
            block: Vec::new(),
            position: 0,
            compressed: 0,
            uncompressed: 0,
            entries: Vec::new(),
        };
        let (index, reader) = fai(format, BufReader::new(reader))?;
        let reader = reader.into_inner();
        if reader.compressed == 0 {
            return Err(invalid("Compressed input is empty, not BGZF"));
        }
        (index, Some(reader.entries))
    } else {
        (fai(format, BufReader::new(file))?.0, None)
    };
    File::create(format!("{path}.fai"))?.write_all(&index)?;
    if let Some(entries) = entries {
        let mut gzi = File::create(format!("{path}.gzi"))?;
        gzi.write_u64::<LittleEndian>(entries.len() as u64)?;
        for (compressed, uncompressed) in entries {
            gzi.write_u64::<LittleEndian>(compressed)?;
            gzi.write_u64::<LittleEndian>(uncompressed)?;
        }
        gzi.flush()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reader_indexed::IndexedReader;
    use bgzip::{BGZFWriter, Compression};
    use std::sync::atomic::{AtomicUsize, Ordering};

    struct Temp(std::path::PathBuf);
    impl Temp {
        fn new() -> Self {
            static ID: AtomicUsize = AtomicUsize::new(0);
            let path = std::env::temp_dir().join(format!(
                "rsbio-index-{}-{}",
                std::process::id(),
                ID.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir_all(&path).unwrap();
            Self(path)
        }
        fn path(&self, name: &str) -> String {
            self.0.join(name).to_str().unwrap().to_string()
        }
    }
    impl Drop for Temp {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn fixture_indexes_match_reference_and_fetch_records() {
        let temp = Temp::new();
        for suffix in ["fa", "fq", "fa.gz", "fq.gz"] {
            let source = format!("test_data/reads_indexable.{suffix}");
            let path = temp.path(&format!("reads.{suffix}"));
            std::fs::copy(&source, &path).unwrap();
            build_index(&path).unwrap();
            assert_eq!(
                std::fs::read(format!("{path}.fai")).unwrap(),
                std::fs::read(format!("{source}.fai")).unwrap()
            );
            let gzi = path.ends_with(".gz").then(|| format!("{path}.gzi"));
            if let Some(gzi) = &gzi {
                assert_eq!(
                    std::fs::read(gzi).unwrap(),
                    std::fs::read(format!("{source}.gzi")).unwrap()
                );
            }
            let mut reader = IndexedReader::new(
                SeqFormat::get(&path).unwrap(),
                &path,
                &format!("{path}.fai"),
                gzi.as_deref(),
            )
            .unwrap();
            let mut sequential = crate::reader::Sequences::new(
                SeqFormat::get(&path).unwrap(),
                crate::reader::get_reader(&path).unwrap(),
            );
            while let Some(record) = sequential.read_next().unwrap() {
                assert_eq!(reader.get_record(&record.id).unwrap(), record);
            }
        }
    }

    #[test]
    fn crlf_wrapping_quality_markers_and_missing_final_newline() {
        let temp = Temp::new();
        for (suffix, input, expected) in [
            (
                "fa",
                b">one\tfirst\r\nAC\r\nGT\r\nN\r\n>two\r\nT".as_slice(),
                "",
            ),
            (
                "fq",
                b"@one\tfirst\r\nAC\r\nGT\r\nN\r\n+one\r\n@+\r\n!I\r\n#\r\n@two\r\nT\r\n+\r\n@"
                    .as_slice(),
                "@+!I#",
            ),
        ] {
            for compressed in [false, true] {
                let path = temp.path(&format!(
                    "reads.{suffix}{}",
                    if compressed { ".gz" } else { "" }
                ));
                if compressed {
                    let mut file = File::create(&path).unwrap();
                    // Force headers and sequence lines to cross block boundaries.
                    for chunk in input.chunks(7) {
                        let mut writer = BGZFWriter::new(&mut file, Compression::default());
                        writer.write_all(chunk).unwrap();
                        writer.close().unwrap();
                    }
                } else {
                    std::fs::write(&path, input).unwrap();
                }
                build_index(&path).unwrap();
                let gzi = compressed.then(|| format!("{path}.gzi"));
                let mut reader = IndexedReader::new(
                    SeqFormat::get(&path).unwrap(),
                    &path,
                    &format!("{path}.fai"),
                    gzi.as_deref(),
                )
                .unwrap();
                let second = reader.get_record("two").unwrap();
                assert_eq!(second.id, "two");
                assert_eq!(second.seq, "T");
                assert_eq!(second.qual, if suffix == "fq" { "@" } else { "" });
                let first = reader.get_record("one").unwrap();
                assert_eq!(first.desc, "first");
                assert_eq!(first.seq, "ACGTN");
                assert_eq!(first.qual, expected);
            }
        }
    }

    #[test]
    fn invalid_input_preserves_existing_indexes() {
        let temp = Temp::new();
        for (suffix, input) in [
            ("fa", b">x\nAA\nA\nAA\n".as_slice()),
            ("fa", b">x\nAA\n>x\nT\n"),
            ("fa", b">empty\n"),
            ("fa", b">\nAA\n"),
            ("fq", b"@x\nAA\n+\n!"),
            ("fq", b"@x\nAA\n+\n!!!\n"),
            ("fq", b"@x\nAA\n"),
            ("fq", b"@x\nAA\nTT\n+\n!!!!\n"),
        ] {
            let path = temp.path(&format!("bad.{suffix}"));
            std::fs::write(&path, input).unwrap();
            std::fs::write(format!("{path}.fai"), b"existing index").unwrap();
            assert!(build_index(&path).is_err());
            assert_eq!(
                std::fs::read(format!("{path}.fai")).unwrap(),
                b"existing index"
            );
        }
    }

    #[test]
    fn ordinary_gzip_corrupt_and_truncated_bgzf_are_rejected() {
        let temp = Temp::new();
        let path = temp.path("bad.fa.gz");
        let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip.write_all(b">x\nAA\n").unwrap();
        std::fs::write(&path, gzip.finish().unwrap()).unwrap();
        assert!(build_index(&path).unwrap_err().to_string().contains("BGZF"));
        let mut data = std::fs::read("test_data/reads_indexable.fa.gz").unwrap();
        data[20] ^= 255;
        std::fs::write(&path, &data).unwrap();
        assert!(build_index(&path).is_err());
        let data = std::fs::read("test_data/reads_indexable.fa.gz").unwrap();
        std::fs::write(&path, &data[..data.len() - 5]).unwrap();
        assert!(build_index(&path).is_err());
        assert!(!std::path::Path::new(&format!("{path}.fai")).exists());
    }
}
