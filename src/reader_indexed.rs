use byteorder::{LittleEndian, ReadBytesExt};
use flate2::read::MultiGzDecoder;

use crate::seq::{SeqFormat, Sequence};
use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, Read, Seek},
};

pub trait SeekRead: Seek + Read + Sync + Send {}

impl SeekRead for File {}

#[derive(Debug, PartialEq)]
pub struct FastaOffsets {
    name: String,
    start: usize,
    end: usize,
    length: usize,
    offset: usize,
    line_bases: usize,
    line_width: usize,
}

pub struct IndexedReader {
    // in future we may want to support fastq
    #[allow(dead_code)]
    format: SeqFormat,
    file: Box<dyn SeekRead>,
    index: HashMap<String, FastaOffsets>,
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
            if parts.len() == 5 {
                let numbers: Result<Vec<usize>, _> =
                    parts[1..].iter().map(|x| x.parse::<usize>()).collect();
                let numbers = numbers.map_err(|_| "Error parsing offset")?;
                if let [length, offset, line_bases, line_width] = numbers[..] {
                    let total_lines = length.div_ceil(line_bases);
                    let total_chars = (total_lines - 1) * line_width
                        + (length - (total_lines - 1) * line_bases)
                        + 1;

                    // we have already checked the length of parts
                    unsafe {
                        index.insert(
                            parts.get_unchecked(0).to_string(),
                            FastaOffsets {
                                start,
                                end: offset + total_chars - 1,
                                name: parts.get_unchecked(0).to_string(),
                                length,
                                offset,
                                line_bases,
                                line_width,
                            },
                        )
                    };
                    start += total_chars + offset;
                } else {
                    return Err(format!("Unrecognised format in line: {}", line));
                }
            } else {
                todo!("Handle other fai formats");
            }
        }

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

        if !self.gzipped {
            let mut seq: Vec<u8> = vec![0; index_entry.length];
            self.file
                .seek(std::io::SeekFrom::Start(index_entry.start as u64))
                .map_err(|_| "Seek error")?;
            let mut id_line = vec![0; index_entry.offset - index_entry.start - 1];
            self.file
                .read_exact(&mut id_line)
                .map_err(|_| "Error reading id line")?;
            let id_line = String::from_utf8_lossy(&id_line).into_owned();

            if let Some((new_id, new_desc)) = id_line.trim_start_matches('>').split_once(' ') {
                id = new_id.to_string();
                desc = new_desc.to_string();
            } else {
                id = id_line;
                desc = String::new();
            }

            self.file
                .seek(std::io::SeekFrom::Start(index_entry.offset as u64))
                .map_err(|_| "Seek error")?;

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

            Ok(Sequence {
                id,
                desc,
                seq: String::from_utf8_lossy(&seq).into_owned(),
                qual: String::new(),
            })
        } else {
            let mut seq = String::new();
            // find start and end of the compressed block
            let gzi_index = self.gzi_index.as_ref().ok_or("No gzi index")?;
            let start_index = gzi_index.binary_search_by(|&(_, uncompressed_offset)| {
                if uncompressed_offset <= index_entry.start as u64 {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Greater
                }
            });
            let start_index = match start_index {
                Ok(start) => start,
                Err(start) => {
                    if start > 0 {
                        start - 1
                    } else {
                        start
                    }
                }
            };
            let end_index = gzi_index.binary_search_by(|&(_, uncompressed_offset)| {
                if uncompressed_offset <= index_entry.end as u64 {
                    std::cmp::Ordering::Greater
                } else {
                    std::cmp::Ordering::Less
                }
            });
            let end_index = match end_index {
                Ok(end) => end,
                Err(end) => {
                    if end < gzi_index.len() - 1 {
                        end + 1
                    } else {
                        gzi_index.len() - 1
                    }
                }
            };

            let mut reader = BufReader::new(self.file.by_ref());
            reader
                .seek(io::SeekFrom::Start(gzi_index[start_index].0))
                .map_err(|_| "Seek error")?;
            let reader = reader.take(gzi_index[end_index].1 - gzi_index[start_index].0);
            let decoder = MultiGzDecoder::new(reader);
            let mut decompressed_reader = BufReader::new(decoder);
            let mut buffer = String::new();
            let mut recording = false;

            while decompressed_reader
                .read_line(&mut buffer)
                .map_err(|_| "Unable to read lines")?
                > 0
            {
                // check if this is the valid fasta header
                if buffer.starts_with('>') {
                    let (new_id, new_desc) = if let Some((new_id, new_desc)) =
                        buffer.trim().trim_start_matches('>').split_once(' ')
                    {
                        (new_id.to_string(), new_desc.to_string())
                    } else {
                        (buffer.trim().to_string(), String::new())
                    };

                    if new_id == name {
                        desc = new_desc;
                        id = new_id;
                        recording = true;
                        buffer.clear();
                        continue;
                    }
                }

                // we have found the record, now we can start recording the sequence
                if recording {
                    // definitely the next entry, break
                    if buffer.starts_with('>') {
                        buffer.clear();
                        break;
                    }
                    seq.push_str(buffer.trim());
                }
                buffer.clear();
            }

            Ok(Sequence {
                id,
                desc,
                seq,
                qual: String::new(),
            })
        }
    }

    pub fn has_record(&self, name: &str) -> bool {
        self.index.contains_key(name)
    }
}

#[cfg(test)]
mod reader_indexed_tests {
    use super::*;
    use crate::seq::SeqFormat;
    const PATH_FA: &str = "test_data/reads_indexable.fa";
    const PATH_FA_GZ: &str = "test_data/reads_indexable.fa.gz";
    const PATH_FA_GZ_IX: &str = "test_data/reads_indexable.fa.gz.fai";
    const PATH_FA_GZ_GZI: &str = "test_data/reads_indexable.fa.gz.gzi";
    const PATH_FA_IX: &str = "test_data/reads_indexable.fa.fai";

    #[test]
    fn fasta_index_load_test() {
        let reader = IndexedReader::new(SeqFormat::Fasta, PATH_FA, PATH_FA_IX, None).unwrap();

        assert_eq!(
            *reader.index.get("Record_1").unwrap(),
            FastaOffsets {
                name: "Record_1".to_string(),
                start: 0,
                end: 313,
                length: 286,
                offset: 24,
                line_bases: 72,
                line_width: 73,
            }
        );
        assert_eq!(
            *reader.index.get("Record_2").unwrap(),
            FastaOffsets {
                name: "Record_2".to_string(),
                start: 314,
                end: 479,
                length: 140,
                offset: 338,
                line_bases: 70,
                line_width: 71,
            }
        );
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
        let entries = reader.gzi_index.unwrap();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries, [(0, 0), (553, 65254)]);
        assert_eq!(
            *reader.index.get("Record_1").unwrap(),
            FastaOffsets {
                name: "Record_1".to_string(),
                start: 0,
                end: 313,
                length: 286,
                offset: 24,
                line_bases: 72,
                line_width: 73,
            }
        );
        assert_eq!(
            *reader.index.get("Record_2").unwrap(),
            FastaOffsets {
                name: "Record_2".to_string(),
                start: 314,
                end: 479,
                length: 140,
                offset: 338,
                line_bases: 70,
                line_width: 71,
            }
        );
    }

    #[test]
    fn fasta_record_load_test() {
        let mut reader = IndexedReader::new(SeqFormat::Fasta, PATH_FA, PATH_FA_IX, None).unwrap();
        let record = reader.get_record("Record_1").unwrap();

        assert_eq!(
            Sequence {
                id: "Record_1".to_string(),
                seq: "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA\
TGCATGCTGATCGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAC\
CGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTG\
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
                    .to_string(),
                desc: "Description 1".to_string(),
                qual: "".to_string(),
            },
            record
        );

        let record = reader.get_record("Record_2").unwrap();

        assert_eq!(
            Sequence {
                id: "Record_2".to_string(),
                seq: "GTTCAGGGATACGACGTTTGTATTTTAAGAATCCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA\
TGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCATCAAGTTACCCTTAACAACTTAAGGGTTTTCA"
                    .to_string(),
                desc: "Description 2".to_string(),
                qual: "".to_string(),
            },
            record
        )
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
        let record = reader.get_record("Record_1").unwrap();

        assert_eq!(
            Sequence {
                id: "Record_1".to_string(),
                seq: "GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA\
TGCATGCTGATCGATCGTACGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAC\
CGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTG\
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
                    .to_string(),
                desc: "Description 1".to_string(),
                qual: "".to_string(),
            },
            record
        );

        let record = reader.get_record("Record_2").unwrap();

        assert_eq!(
            Sequence {
                id: "Record_2".to_string(),
                seq: "GTTCAGGGATACGACGTTTGTATTTTAAGAATCCAAGTTACCCTTAACAACTTAAGGGTTTTCAAATAGA\
        TGAAGCAGAAGTCGATGATAATACGCGTCGTTTTATCATCAAGTTACCCTTAACAACTTAAGGGTTTTCA"
                    .to_string(),
                desc: "Description 2".to_string(),
                qual: "".to_string(),
            },
            record
        )
    }
}
