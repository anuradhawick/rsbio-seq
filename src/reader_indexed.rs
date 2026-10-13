use crate::seq::{SeqFormat, Sequence};
use byteorder::{LittleEndian, ReadBytesExt};
use flate2::read::MultiGzDecoder;
use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, Read, Seek},
};

#[derive(Debug, PartialEq)]
pub struct FastXOffsets {
    name: String,
    start: usize,
    end: usize,
    length: usize,
    offset: usize,
    line_bases: usize,
    line_width: usize,
    qual_offset: Option<usize>,
}

pub struct IndexedReader {
    file: File,
    index: HashMap<String, FastXOffsets>,
    gzi_index: Option<Vec<(u64, u64)>>,
}

impl IndexedReader {
    pub fn new(
        format: SeqFormat,
        records_file: &str,
        index_file: &str,
        gzi_file: Option<&str>,
    ) -> Result<Self, String> {
        let file = File::open(records_file).map_err(|e| e.to_string())?;
        let reader = BufReader::new(File::open(index_file).map_err(|e| e.to_string())?);
        let mut index = HashMap::new();
        let mut start = 0;
        for line in reader.lines() {
            let line = line.map_err(|e| e.to_string())?;
            let parts: Vec<_> = line.split('\t').collect();
            let fields = match format {
                SeqFormat::Fasta => 5,
                SeqFormat::Fastq => 6,
            };
            if parts.len() != fields || parts[0].is_empty() {
                return Err("Invalid FAI index entry".into());
            }
            let numbers = parts[1..]
                .iter()
                .map(|x| x.parse::<usize>())
                .collect::<Result<Vec<_>, _>>()
                .map_err(|_| "Invalid FAI offset")?;
            let (length, offset, line_bases, line_width) =
                (numbers[0], numbers[1], numbers[2], numbers[3]);
            if length == 0 || line_bases == 0 || line_width < line_bases || offset <= start {
                return Err("Invalid FAI line geometry or offset".into());
            }
            let span = ((length - 1) / line_bases)
                .checked_mul(line_width)
                .and_then(|n| n.checked_add((length - 1) % line_bases + 1))
                .ok_or("FAI offset overflow")?;
            let qual_offset = numbers.get(4).copied();
            if qual_offset.is_some_and(|q| q <= offset) {
                return Err("Invalid FAI quality offset".into());
            }
            let end = qual_offset
                .unwrap_or(offset)
                .checked_add(span)
                .ok_or("FAI offset overflow")?;
            let next_start = end
                .checked_add(line_width - line_bases)
                .ok_or("FAI offset overflow")?;
            let name = parts[0].to_string();
            if index
                .insert(
                    name.clone(),
                    FastXOffsets {
                        name,
                        start,
                        end,
                        length,
                        offset,
                        line_bases,
                        line_width,
                        qual_offset,
                    },
                )
                .is_some()
            {
                return Err("Duplicate identifier in FAI index".into());
            }
            start = next_start;
        }
        let gzi_index = if records_file.ends_with(".gz") {
            let mut probe = file.try_clone().map_err(|e| e.to_string())?;
            let mut header = [0; 12];
            probe.read_exact(&mut header).map_err(|e| e.to_string())?;
            if header[..4] != [31, 139, 8, 4] {
                return Err("Gzipped file is not BGZF".into());
            }
            let mut extra = vec![0; u16::from_le_bytes([header[10], header[11]]) as usize];
            probe.read_exact(&mut extra).map_err(|e| e.to_string())?;
            let mut p = 0;
            let mut bgzf = false;
            while p + 4 <= extra.len() {
                let n = u16::from_le_bytes([extra[p + 2], extra[p + 3]]) as usize;
                if p + 4 + n > extra.len() {
                    return Err("Invalid BGZF header".into());
                }
                bgzf |= &extra[p..p + 2] == b"BC" && n == 2;
                p += 4 + n;
            }
            if !bgzf || p != extra.len() {
                return Err("Gzipped file is not BGZF".into());
            }
            let mut reader = BufReader::new(
                File::open(gzi_file.ok_or("BGZF input requires a GZI index")?)
                    .map_err(|e| e.to_string())?,
            );
            let count = reader
                .read_u64::<LittleEndian>()
                .map_err(|e| e.to_string())?;
            let mut entries = vec![(0, 0)];
            for _ in 0..count {
                let compressed = reader
                    .read_u64::<LittleEndian>()
                    .map_err(|e| e.to_string())?;
                let uncompressed = reader
                    .read_u64::<LittleEndian>()
                    .map_err(|e| e.to_string())?;
                let previous = entries.last().unwrap();
                if compressed <= previous.0 || uncompressed <= previous.1 {
                    return Err("GZI offsets must be increasing".into());
                }
                entries.push((compressed, uncompressed));
            }
            Some(entries)
        } else {
            None
        };
        Ok(Self {
            file,
            index,
            gzi_index,
        })
    }

    // FAI offsets refer to uncompressed bytes. Start decoding at the nearest
    // preceding BGZF block, then skip only within that block.
    fn read_region(&mut self, offset: usize, length: usize) -> Result<Vec<u8>, String> {
        let mut bytes = vec![0; length];
        if let Some(entries) = &self.gzi_index {
            let i = entries
                .partition_point(|&(_, u)| u <= offset as u64)
                .saturating_sub(1);
            let (compressed, uncompressed) = entries[i];
            self.file
                .seek(io::SeekFrom::Start(compressed))
                .map_err(|e| e.to_string())?;
            let mut decoder = MultiGzDecoder::new(&mut self.file);
            let skip = offset as u64 - uncompressed;
            let copied = io::copy(&mut decoder.by_ref().take(skip), &mut io::sink())
                .map_err(|e| e.to_string())?;
            if copied != skip {
                return Err("Truncated BGZF input".into());
            }
            decoder.read_exact(&mut bytes).map_err(|e| e.to_string())?;
        } else {
            self.file
                .seek(io::SeekFrom::Start(offset as u64))
                .map_err(|e| e.to_string())?;
            self.file
                .read_exact(&mut bytes)
                .map_err(|e| e.to_string())?;
        }
        Ok(bytes)
    }

    pub fn get_record(&mut self, name: &str) -> Result<Sequence, String> {
        let entry = self.index.get(name).ok_or("Record not found")?;
        let (start, offset, length, bases, width, qual_offset) = (
            entry.start,
            entry.offset,
            entry.length,
            entry.line_bases,
            entry.line_width,
            entry.qual_offset,
        );
        let header = self.read_region(start, offset - start)?;
        let header = std::str::from_utf8(&header)
            .map_err(|e| e.to_string())?
            .trim_end_matches(['\r', '\n']);
        if !header.starts_with(['>', '@']) {
            return Err("Invalid record header".into());
        }
        let header = &header[1..];
        let split = header.find(char::is_whitespace).unwrap_or(header.len());
        let id = header[..split].to_string();
        let desc = if split < header.len() {
            header[split + 1..].to_string()
        } else {
            String::new()
        };
        if id != name {
            return Err("FAI index does not match record header".into());
        }
        let span = ((length - 1) / bases) * width + (length - 1) % bases + 1;
        let seq = self.read_bases(offset, span, length)?;
        let qual = if let Some(offset) = qual_offset {
            self.read_bases(offset, span, length)?
        } else {
            String::new()
        };
        Ok(Sequence {
            id,
            desc,
            seq,
            qual,
        })
    }

    fn read_bases(&mut self, offset: usize, span: usize, length: usize) -> Result<String, String> {
        let mut bytes = self.read_region(offset, span)?;
        bytes.retain(|&b| b != b'\r' && b != b'\n');
        if bytes.len() != length {
            return Err("FAI index does not match sequence line geometry".into());
        }
        String::from_utf8(bytes).map_err(|e| e.to_string())
    }

    pub fn has_record(&self, name: &str) -> bool {
        self.index.contains_key(name)
    }
}

#[cfg(test)]
mod reader_indexed_fasta_tests {
    use super::*;
    use crate::seq::SeqFormat;
    const PATH_FA: &str = "test_data/reads_indexable.fa";
    const PATH_FA_GZ: &str = "test_data/reads_indexable.fa.gz";
    const PATH_FA_GZ_IX: &str = "test_data/reads_indexable.fa.gz.fai";
    const PATH_FA_GZ_GZI: &str = "test_data/reads_indexable.fa.gz.gzi";
    const PATH_FA_IX: &str = "test_data/reads_indexable.fa.fai";

    fn assert_index_entries(reader: &IndexedReader) {
        assert_eq!(reader.index.len(), 7);
        assert_eq!(
            *reader.index.get("Record_1").unwrap(),
            FastXOffsets {
                name: "Record_1".to_string(),
                start: 0,
                end: 313,
                length: 286,
                offset: 24,
                line_bases: 72,
                line_width: 73,
                qual_offset: None,
            }
        );
        assert_eq!(
            *reader.index.get("Record_2").unwrap(),
            FastXOffsets {
                name: "Record_2".to_string(),
                start: 314,
                end: 479,
                length: 140,
                offset: 338,
                line_bases: 70,
                line_width: 71,
                qual_offset: None,
            }
        );
    }

    fn assert_records(reader: &mut IndexedReader) {
        // Use the sequential FASTA parser as an independent sequence reference.
        let records: Vec<_> = crate::reader::Sequences::new(
            SeqFormat::Fasta,
            crate::reader::get_reader("test_data/reads_indexable.fa").unwrap(),
        )
        .collect();

        // Read out of order, including records crossing BGZF block boundaries
        // and repeated lookups to check that each read seeks independently.
        for index in [0, 1, 6, 2, 5, 4, 3, 0] {
            let expected = &records[index];
            assert_eq!(
                reader.get_record(expected.id.as_str()).unwrap(),
                Sequence {
                    id: expected.id.as_str().to_string(),
                    desc: expected.desc.as_str().to_string(),
                    seq: expected.seq.clone(),
                    qual: String::new(),
                }
            );
        }
        assert!(reader.has_record("Record_1"));
        assert!(!reader.has_record("DoesNotExist"));
        assert_eq!(
            reader.get_record("DoesNotExist"),
            Err("Record not found".to_string())
        );
        assert_eq!(reader.get_record("Record_2").unwrap().id, "Record_2");
    }

    #[test]
    fn fasta_index_load_test() {
        let reader = IndexedReader::new(SeqFormat::Fasta, PATH_FA, PATH_FA_IX, None).unwrap();
        assert_index_entries(&reader);
        assert!(reader.gzi_index.is_none());
    }

    #[test]
    fn fasta_index_load_gz_test() {
        let reader = IndexedReader::new(
            SeqFormat::Fasta,
            PATH_FA_GZ,
            PATH_FA_GZ_IX,
            Some(PATH_FA_GZ_GZI),
        )
        .unwrap();
        assert_index_entries(&reader);
        assert_eq!(
            reader.gzi_index.as_ref().unwrap(),
            &vec![(0, 0), (553, 65254)]
        );
    }

    #[test]
    fn fasta_record_load_test() {
        let mut reader = IndexedReader::new(SeqFormat::Fasta, PATH_FA, PATH_FA_IX, None).unwrap();
        assert_records(&mut reader);
    }

    #[test]
    fn fasta_record_load_gz_test() {
        let mut reader = IndexedReader::new(
            SeqFormat::Fasta,
            PATH_FA_GZ,
            PATH_FA_GZ_IX,
            Some(PATH_FA_GZ_GZI),
        )
        .unwrap();
        assert_records(&mut reader);
    }
}

#[cfg(test)]
mod reader_indexed_fastq_tests {
    use super::*;
    use crate::seq::SeqFormat;
    const PATH_FQ: &str = "test_data/reads_indexable.fq";
    const PATH_FQ_GZ: &str = "test_data/reads_indexable.fq.gz";
    const PATH_FQ_GZ_IX: &str = "test_data/reads_indexable.fq.gz.fai";
    const PATH_FQ_GZ_GZI: &str = "test_data/reads_indexable.fq.gz.gzi";
    const PATH_FQ_IX: &str = "test_data/reads_indexable.fq.fai";

    fn assert_index_entries(reader: &IndexedReader) {
        assert_eq!(reader.index.len(), 7);
        assert_eq!(
            *reader.index.get("Record_1").unwrap(),
            FastXOffsets {
                name: "Record_1".to_string(),
                start: 0,
                end: 599,
                length: 286,
                offset: 24,
                line_bases: 286,
                line_width: 287,
                qual_offset: Some(313),
            }
        );
        assert_eq!(
            *reader.index.get("Record_2").unwrap(),
            FastXOffsets {
                name: "Record_2".to_string(),
                start: 600,
                end: 907,
                length: 140,
                offset: 624,
                line_bases: 140,
                line_width: 141,
                qual_offset: Some(767),
            }
        );
    }

    fn assert_records(reader: &mut IndexedReader) {
        // Use the matching FASTA fixture as an independent sequence reference.
        let records: Vec<_> = crate::reader::Sequences::new(
            SeqFormat::Fasta,
            crate::reader::get_reader("test_data/reads_indexable.fa").unwrap(),
        )
        .collect();

        // Read out of order, including records crossing BGZF block boundaries
        // and repeated lookups to check that each read seeks independently.
        for index in [0, 1, 6, 2, 5, 4, 3, 0] {
            let expected = &records[index];
            let length = expected.seq.len();
            let qual = "!5I@+".repeat(length.div_ceil(5));
            assert_eq!(
                reader.get_record(expected.id.as_str()).unwrap(),
                Sequence {
                    id: expected.id.as_str().to_string(),
                    desc: expected.desc.as_str().to_string(),
                    seq: expected.seq.clone(),
                    qual: qual[..length].to_string(),
                }
            );
        }
        assert!(reader.has_record("Record_1"));
        assert!(!reader.has_record("DoesNotExist"));
        assert_eq!(
            reader.get_record("DoesNotExist"),
            Err("Record not found".to_string())
        );
        assert_eq!(reader.get_record("Record_2").unwrap().id, "Record_2");
    }

    #[test]
    fn fastq_index_load_test() {
        let reader = IndexedReader::new(SeqFormat::Fastq, PATH_FQ, PATH_FQ_IX, None).unwrap();
        assert_index_entries(&reader);
        assert!(reader.gzi_index.is_none());
    }

    #[test]
    fn fastq_index_load_gz_test() {
        let reader = IndexedReader::new(
            SeqFormat::Fastq,
            PATH_FQ_GZ,
            PATH_FQ_GZ_IX,
            Some(PATH_FQ_GZ_GZI),
        )
        .unwrap();
        assert_index_entries(&reader);
        assert_eq!(
            reader.gzi_index.as_ref().unwrap(),
            &vec![(0, 0), (55, 24), (449, 46295), (718, 105991)]
        );
    }

    #[test]
    fn fastq_record_load_test() {
        let mut reader = IndexedReader::new(SeqFormat::Fastq, PATH_FQ, PATH_FQ_IX, None).unwrap();
        assert_records(&mut reader);
    }

    #[test]
    fn fastq_record_load_gz_test() {
        let mut reader = IndexedReader::new(
            SeqFormat::Fastq,
            PATH_FQ_GZ,
            PATH_FQ_GZ_IX,
            Some(PATH_FQ_GZ_GZI),
        )
        .unwrap();
        assert_records(&mut reader);
    }
}
