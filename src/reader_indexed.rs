use byteorder::{LittleEndian, ReadBytesExt};
use flate2::bufread::GzDecoder;

use crate::seq::{SeqFormat, Sequence};
use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, Read, Seek},
};

pub trait SeekRead: Seek + Read + Sync + Send {}

impl SeekRead for File {}

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
    // in future we may want to support fastq
    #[allow(dead_code)]
    format: SeqFormat,
    file: Box<dyn SeekRead>,
    index: HashMap<String, FastXOffsets>,
    gzipped: bool,
    gzi_index: Option<Vec<(u64, u64)>>,
}

fn is_bgzf(file_path: &str) -> io::Result<bool> {
    let mut file = File::open(file_path)?;
    let mut header = [0u8; 18]; // Read the first 18 bytes of the header
    file.read_exact(&mut header)?;

    // Check for GZIP magic number and compression method
    if header[0] == 0x1f && header[1] == 0x8b && header[2] == 0x08 {
        // Check for BGZF-specific bytes at position 12 and 13
        if header[12] == 0x42 && header[13] == 0x43 {
            return Ok(true);
        }
    }
    Ok(false)
}

impl IndexedReader {
    pub fn new(
        format: SeqFormat,
        records_file: &str,
        index_file: &str,
        gzi_file: Option<&str>,
    ) -> Result<Self, String> {
        let gzipped = records_file.ends_with(".gz");
        let mut gzi_index = None;
        let file =
            File::open(records_file).map_err(|_| format!("Unable to open: {}", records_file))?;
        let mut index = HashMap::new();
        let index_file =
            File::open(index_file).map_err(|_| format!("Unable to open: {}", index_file))?;
        let reader = std::io::BufReader::new(index_file);
        let mut start = 0;
        for line in reader.lines() {
            let line = line.map_err(|_| "Error reading index file")?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 5 {
                let numbers: Result<Vec<usize>, _> =
                    parts[1..].iter().map(|x| x.parse::<usize>()).collect();
                let numbers = numbers.map_err(|_| "Error parsing offset")?;
                if let [length, offset, line_bases, line_width] = numbers[..4] {
                    let total_lines = length.div_ceil(line_bases);
                    let total_chars = (total_lines - 1) * line_width
                        + (length - (total_lines - 1) * line_bases)
                        + 1;
                    let end = if let Some(qual_offset) = numbers.get(4).copied() {
                        qual_offset + total_chars - 1
                    } else {
                        offset + total_chars - 1
                    };

                    // we have already checked the length of parts
                    unsafe {
                        index.insert(
                            parts.get_unchecked(0).to_string(),
                            FastXOffsets {
                                start,
                                end,
                                name: parts.get_unchecked(0).to_string(),
                                length,
                                offset,
                                line_bases,
                                line_width,
                                qual_offset: numbers.get(4).copied(),
                            },
                        )
                    };
                    println!(
                        "Total chars {}, start {}, offset {}",
                        total_chars, start, offset
                    );
                    if let Some(qual_offset) = numbers.get(4).copied() {
                        start = qual_offset + total_chars;
                    } else {
                        start = offset + total_chars;
                    }
                } else {
                    return Err(format!("Unrecognised format in line: {}", line));
                }
            } else {
                todo!("Handle other fai formats");
            }
        }

        println!("Index: {:#?}", index);

        if gzipped {
            if let Ok(bgzf) = is_bgzf(records_file) {
                if !bgzf {
                    return Err("Gzipped file is not BGZF".to_string());
                }
            }
            let gzi_file = gzi_file.ok_or("Gzipped file requires index")?;
            let gzi_file = File::open(gzi_file).map_err(|_| "Unable to open gzi file")?;
            let mut reader = BufReader::new(gzi_file);
            let no_entries = reader
                .read_u64::<LittleEndian>()
                .map_err(|_| "Error reading gzi")?;
            let mut entries: Vec<(u64, u64)> = Vec::with_capacity(no_entries as usize + 1);
            entries.push((0, 0));
            for _ in 0..no_entries {
                let start = reader
                    .read_u64::<LittleEndian>()
                    .map_err(|_| "Error reading gzi")?;
                let end = reader
                    .read_u64::<LittleEndian>()
                    .map_err(|_| "Error reading gzi")?;
                entries.push((start, end));
            }

            gzi_index = Some(entries);
        }

        Ok(Self {
            format,
            gzipped,
            file: Box::new(file),
            index,
            gzi_index,
        })
    }

    pub fn get_record(&mut self, name: &str) -> Result<Sequence, String> {
        let index_entry = self.index.get(name).ok_or("Record not found")?;
        let mut buffer: Vec<u8> = vec![0; index_entry.line_width];
        let mut bytes_read = 0;
        let mut id: String = String::new();
        let mut desc: String = String::new();

        println!("Fetching record: {}", name);
        println!("Index entry: {:?}", index_entry);

        if !self.gzipped {
            let mut seq: Vec<u8> = vec![0; index_entry.length];
            let mut qual: Vec<u8> = vec![
                0;
                if index_entry.qual_offset.is_some() {
                    index_entry.length
                } else {
                    0
                }
            ];
            self.file
                .seek(std::io::SeekFrom::Start(index_entry.start as u64))
                .map_err(|_| "Seek error")?;
            let mut id_line = vec![0; index_entry.offset - index_entry.start - 1];
            self.file
                .read_exact(&mut id_line)
                .map_err(|_| "Error reading id line")?;
            let id_line = String::from_utf8_lossy(&id_line).into_owned();

            // Parse the ID and description from the header line
            if id_line.starts_with('>') {
                if let Some((new_id, new_desc)) = id_line.trim_start_matches('>').split_once(' ') {
                    id = new_id.to_string();
                    desc = new_desc.to_string();
                } else {
                    id = id_line;
                    desc = String::new();
                }
            } else if id_line.starts_with('@') {
                // Fastq header line starts with '@'
                if let Some((new_id, new_desc)) = id_line.trim_start_matches('@').split_once(' ') {
                    id = new_id.to_string();
                    desc = new_desc.to_string();
                } else {
                    id = id_line;
                    desc = String::new();
                }
            } else {
                return Err("Invalid header line".into());
            }

            // Seek to the start of the sequence data
            self.file
                .seek(std::io::SeekFrom::Start(index_entry.offset as u64))
                .map_err(|_| "Seek error")?;

            // Read the sequence data in chunks according to the line width and line bases
            while bytes_read < index_entry.length {
                // bytes that must be read in the line
                let to_read =
                    std::cmp::min(index_entry.line_width, index_entry.length - bytes_read);
                let mut handle = self.file.by_ref().take(to_read as u64);

                handle
                    .read(&mut buffer)
                    .map_err(|_| "Error reading chunk")?;
                // bases that must be copied to the sequence
                let to_copy =
                    std::cmp::min(index_entry.line_bases, index_entry.length - bytes_read);
                seq[bytes_read..bytes_read + to_copy].copy_from_slice(&buffer[..to_copy]);
                bytes_read += to_copy;
            }

            // If the sequence has a quality score, read it as well
            if let Some(qual_offset) = index_entry.qual_offset {
                self.file
                    .seek(std::io::SeekFrom::Start(qual_offset as u64))
                    .map_err(|_| "Seek error")?;

                bytes_read = 0;
                while bytes_read < index_entry.length {
                    // bytes that must be read in the line
                    let to_read =
                        std::cmp::min(index_entry.line_width, index_entry.length - bytes_read);
                    let mut handle = self.file.by_ref().take(to_read as u64);

                    handle
                        .read(&mut buffer)
                        .map_err(|_| "Error reading chunk")?;
                    // bases that must be copied to the qual
                    let to_copy =
                        std::cmp::min(index_entry.line_bases, index_entry.length - bytes_read);
                    qual[bytes_read..bytes_read + to_copy].copy_from_slice(&buffer[..to_copy]);
                    bytes_read += to_copy;
                }
            }

            Ok(Sequence {
                id,
                desc,
                seq: String::from_utf8_lossy(&seq).into_owned(),
                qual: String::from_utf8_lossy(&qual).into_owned(),
            })
        } else {
            let mut seq = String::new();
            let mut qual = String::new();
            // find start and end of the compressed block
            let gzi_index = self.gzi_index.as_ref().ok_or("No gzi index")?;
            let start_index = gzi_index
                .partition_point(|&(_, u)| u <= index_entry.start as u64)
                .saturating_sub(1);
            let end_index = (gzi_index.partition_point(|&(_, u)| u <= index_entry.end as u64) + 1)
                .min(gzi_index.len() - 1);

            // Create reader and move to offset in compressed file
            let mut reader = BufReader::new(self.file.by_ref());
            reader
                .seek(io::SeekFrom::Start(gzi_index[start_index].0))
                .map_err(|_| "Seek error")?;
            let mut decompressed_buffer = String::new();

            // Decompress each block in the range and append to the decompressed buffer
            for (comp, uncomp) in &gzi_index[start_index..=end_index] {
                println!("Block: {:?} {:?}", comp, uncomp);
                let mut block_decoder = GzDecoder::new(&mut reader);
                let mut block_decompressed_buffer = String::new();
                block_decoder
                    .read_to_string(&mut block_decompressed_buffer)
                    .map_err(|_| "Block decompression error")?;
                decompressed_buffer.push_str(&block_decompressed_buffer);
            }

            // Calculate the slice of the decompressed buffer that corresponds to the requested sequence
            let decompressed_offset = index_entry.start as u64 - gzi_index[start_index].1;
            let decompressed_end =
                decompressed_offset + (index_entry.end as u64 - index_entry.start as u64);
            let decompressed_slice =
                &decompressed_buffer[decompressed_offset as usize..decompressed_end as usize];

            println!("Decompressed slice: {:?}", decompressed_slice);

            let mut buffer = String::new();
            let mut recording = false;
            let mut recording_qual = false;

            for line in decompressed_slice.lines() {
                // check if this is the valid fasta header
                if line.starts_with('>') {
                    let (new_id, new_desc) = if let Some((new_id, new_desc)) =
                        line.trim().trim_start_matches('>').split_once(' ')
                    {
                        (new_id.to_string(), new_desc.to_string())
                    } else {
                        (line.trim().to_string(), String::new())
                    };

                    if new_id == name {
                        desc = new_desc;
                        id = new_id;
                        recording = true;
                        recording_qual = false;
                        buffer.clear();
                        continue;
                    }
                } else if line.starts_with('@') {
                    let (new_id, new_desc) = if let Some((new_id, new_desc)) =
                        line.trim().trim_start_matches('@').split_once(' ')
                    {
                        (new_id.to_string(), new_desc.to_string())
                    } else {
                        (line.trim().to_string(), String::new())
                    };

                    if new_id == name {
                        desc = new_desc;
                        id = new_id;
                        recording = true;
                        recording_qual = false;
                        buffer.clear();
                        continue;
                    }
                } else if line.starts_with('+') && line.trim().len() == 1 {
                    recording_qual = true;
                    recording = false;
                    buffer.clear();
                    continue;
                }

                // we have found the record, now we can start recording the sequence
                if recording {
                    // definitely the next entry, break (only possible for fasta)
                    if line.starts_with('>') {
                        buffer.clear();
                        break;
                    }
                    seq.push_str(line.trim());
                }
                if recording_qual {
                    // definitely the next entry, break (only possible for fastq)
                    if line.starts_with('@') {
                        buffer.clear();
                        break;
                    }
                    qual.push_str(line.trim());
                }
            }

            Ok(Sequence {
                id,
                desc,
                seq,
                qual,
            })
        }
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
