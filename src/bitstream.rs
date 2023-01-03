use log::debug;
use std::cell::RefCell;
use std::fmt;
use std::path::Path;

use crate::decoder::CodecId;
use num_enum::{FromPrimitive, IntoPrimitive};

pub mod reader;

#[derive(Copy, Clone, Default, Debug)]
pub struct Position {
    pub bytes: usize,
    pub bits: u8,
}

pub struct GofStat {
    pub v3c_unit_size: Vec<usize>,
    pub video_bin_size: Vec<usize>,
}

pub struct Stat {
    pub header: usize,
    pub bitstream_gof_stat: Vec<GofStat>,
}

impl Stat {
    pub fn new() -> Self {
        Self {
            header: 0,
            bitstream_gof_stat: Vec::new(),
        }
    }

    pub fn incr_header(&mut self, size: usize) {
        self.header += size;
    }

    pub fn new_gof(&mut self) {
        self.bitstream_gof_stat.push(GofStat {
            v3c_unit_size: Vec::new(),
            video_bin_size: Vec::new(),
        });
    }
}

#[derive(Default, Clone)]
pub struct Bitstream {
    pub data: Vec<u8>,
    pub position: RefCell<Position>,
}

impl Bitstream {
    pub fn new(data: Vec<u8>) -> Self {
        Self {
            data,
            position: RefCell::new(Position { bytes: 0, bits: 0 }),
        }
    }

    pub fn from_file(path: &Path) -> Self {
        let data = std::fs::read(path).unwrap();
        Self::new(data)
    }

    /// Clear all data in bitstream
    pub fn clear(&mut self) {
        self.data.clear();
        self.reset();
    }

    /// Resets the position to the beginning of the bitstream
    pub fn reset(&self) {
        self.position.borrow_mut().bytes = 0;
        self.position.borrow_mut().bits = 0;
    }

    fn bytes(&self) -> usize {
        self.position.borrow().bytes
    }

    fn bits(&self) -> u8 {
        self.position.borrow().bits
    }

    /// Originally: Bitstream::size
    fn position_in_bytes(&self) -> usize {
        self.position.borrow().bytes
    }

    /// Originally: Bitstream::capacity
    pub fn size(&self) -> u64 {
        self.data.len() as u64
    }

    pub fn more_data(&self) -> bool {
        self.position_in_bytes() < self.data.len()
    }

    #[inline]
    fn is_byte_aligned(&self) -> bool {
        self.bits() == 0
    }

    #[inline]
    fn move_to_next_byte(&self) {
        self.position.borrow_mut().bytes += 1;
        self.position.borrow_mut().bits = 0;
    }

    /// 8.3.3. Byte alignment syntax <==> 8.3.6.10 RBSP trailing bit syntax
    pub fn byte_align(&self) {
        // FIXME: This is a strange wrinkle from the tmc2 repo have the following line
        self.read(1);
        if !self.is_byte_aligned() {
            self.position.borrow_mut().bits = 0;
            self.position.borrow_mut().bytes += 1;
        }
    }

    /// Copy bitstream from src to self
    /// At the end of copy_from, the position of self is updated to the end of the copied data
    pub fn copy_from(&mut self, src: &Bitstream, start_byte: usize, bitstream_size: usize) {
        if self.data.len() < self.bytes() + bitstream_size {
            self.data.resize(self.bytes() + bitstream_size, 0)
        }

        let bytes = self.bytes();
        self.data[bytes..bytes + bitstream_size]
            .copy_from_slice(&src.data[start_byte..start_byte + bitstream_size]);

        self.position.borrow_mut().bytes += bitstream_size;
        src.position.borrow_mut().bytes += bitstream_size;
    }

    pub fn read(&self, bits: u8) -> u32 {
        if bits > 32 {
            panic!("Bitstream::read: bits > 32");
        }

        let mut val: u32 = 0;
        let mut i = 0;
        while i < bits {
            val |= (((self.data[self.bytes()] >> (7 - self.bits())) & 1) as u32) << (bits - i - 1);
            self.position.borrow_mut().bits += 1;
            if self.bits() == 8 {
                self.position.borrow_mut().bytes += 1;
                self.position.borrow_mut().bits = 0;
            }
            i += 1;
        }
        val
    }

    pub fn peek(&self, bits: u8) -> u32 {
        let pos = self.position.borrow().clone();
        let res = self.read(bits);
        // rewind the cursor
        *self.position.borrow_mut() = pos;
        res
    }

    fn read_slice(&self, size: usize) -> &[u8] {
        let bytes = self.bytes();
        self.position.borrow_mut().bytes += size;
        &self.data[bytes..bytes + size]
    }

    /// read unsigned variable length code (Oth order Exp-Golomb codes)
    fn read_uvlc(&self) -> u32 {
        let mut leading_zeros = 0;
        while self.read(1) == 0 {
            leading_zeros += 1;
        }
        if leading_zeros == 0 {
            return 0;
        }
        (1 << leading_zeros) - 1 + self.read(leading_zeros)
    }

    fn read_svlc(&self) -> i32 {
        let x = self.read_uvlc();
        if x & 1 == 1 {
            (x >> 1) as i32 + 1
        } else {
            -1 * (x >> 1) as i32
        }
    }
}

pub(crate) struct VideoBitstream {
    pub(crate) data: Vec<u8>,
    pub(crate) video_type: VideoType,
}

impl VideoBitstream {
    pub fn from_bitstream(bitstream: &Bitstream, size: usize, video_type: VideoType) -> Self {
        let mut vbs = Self {
            data: vec![0; size],
            video_type,
        };

        vbs.data[..size].copy_from_slice(bitstream.read_slice(size));
        debug!("VideoType={:?} size={}", &vbs.video_type, &vbs.data.len());
        vbs
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Sample stream is used in the final bitstream. However, the decoder only understands the bytestream
    /// so this fn is used in PCCDecoder before being send to the `VideoDecoder::decode`
    pub fn sample_stream_to_bytestream(&self, codec_id: CodecId, precision: usize) -> Vec<u8> {
        // let emulation_prevention_bytes = false;
        // let change_start_code_size = true;
        assert_eq!(precision, 4);

        let mut size_start_code = 4;
        let mut start_index = 0;
        let mut end_index = 0;
        let mut new_frame = true;
        let mut result = Vec::with_capacity(self.data.len());

        loop {
            let mut nalu_size = 0;
            for i in 0..precision {
                nalu_size = (nalu_size << 8) + self.data[start_index + i] as usize;
            }

            end_index = start_index + precision + nalu_size;
            if end_index >= self.data.len() {
                break;
            }

            for i in 0..size_start_code - 1 {
                result.push(0);
            }
            result.push(1);

            // if emulation_prevention_bytes {
            //     unimplemented!("emulation_prevention_bytes");
            // } else {
            for i in (start_index + precision)..end_index {
                result.push(self.data[i]);
            }
            // }

            start_index = end_index;
            if start_index + precision < self.data.len() {
                let mut nalu_type = 0;
                let mut use_long_start_code = false;
                new_frame = false;
                match codec_id {
                    CodecId::H264 => use_long_start_code = true,
                    CodecId::H265 => {
                        nalu_type = (self.data[start_index + precision + 1] & 126) >> 1;
                        use_long_start_code = new_frame || (nalu_type >= 32 && nalu_type < 41);
                        if nalu_type < 12 {
                            new_frame = true;
                        }
                    }
                    CodecId::H266 => {
                        nalu_type = (self.data[start_index + precision + 1] & 248) >> 3;
                        use_long_start_code = new_frame || (nalu_type >= 12 && nalu_type < 20);
                        if nalu_type < 12 {
                            new_frame = true;
                        }
                    }
                }
                size_start_code = if use_long_start_code { 4 } else { 3 };
            }
        }

        result
    }
}

// TODO: set Attribute enum values
#[derive(Debug, PartialEq, FromPrimitive, IntoPrimitive)]
#[repr(u8)]
pub(crate) enum VideoType {
    #[default]
    Occupancy,
    Geometry,
    GeometryD0,
    GeometryD1,
    GeometryD2,
    GeometryD3,
    GeometryD4,
    GeometryD5,
    GeometryD6,
    GeometryD7,
    GeometryD8,
    GeometryD9,
    GeometryD10,
    GeometryD11,
    GeometryD12,
    GeometryD13,
    GeometryD14,
    GeometryD15,
    GeometryRaw,
    Attribute,
    // AttributeT0,
    // AttributeT1 = 20 + 64,
    // AttributeT2,
    // AttributeT3,
    // AttributeT4,
    // AttributeT5,
    // AttributeT6,
    // AttributeT7,
    // AttributeT8,
    // AttributeT9,
    // AttributeT10,
    // AttributeT11,
    // AttributeT12,
    // AttributeT13,
    // AttributeT14,
    // AttributeT15,
    // AttributeRaw,
    // NumVideoType,
}

impl fmt::Display for VideoType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
        // or, alternatively:
        // fmt::Debug::fmt(self, f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitstream_read() {
        let bitstream = Bitstream::new(vec![0b10101010, 0b11110000, 0b11001001, 0b00110011]);
        assert_eq!(bitstream.read(1), 0b1);
        assert_eq!(bitstream.read(3), 0b010);
        assert_eq!(bitstream.read(7), 0b1010111);
        assert_eq!(bitstream.read(11), 0b10000110010);
        assert_eq!(bitstream.read(4), 0b0100);
        assert_eq!(bitstream.read(6), 0b110011);
        bitstream.reset();
        assert_eq!(bitstream.read(8), 0b10101010);
    }

    #[test]
    fn test_bitstream_peek() {
        let bitstream = Bitstream::new(vec![0b10101010]);
        assert_eq!(bitstream.peek(1), 0b1);
        assert_eq!(bitstream.peek(1), 0b1);
        assert_eq!(bitstream.peek(3), 0b101);
        assert_eq!(bitstream.peek(3), 0b101);
    }

    #[test]
    fn test_bitstream_read_uvlc() {
        let bitstream = Bitstream::new(vec![
            0b10100110, 0b01000010, 0b10011000, 0b11100010, 0b00000100, 0b10001010, 0b00010110,
            0b00110000, 0b01101000, 0b11100001, 0b11100000,
        ]);
        assert_eq!(bitstream.read_uvlc(), 0);
        assert_eq!(bitstream.read_uvlc(), 1);
        assert_eq!(bitstream.read_uvlc(), 2);
        assert_eq!(bitstream.read_uvlc(), 3);
        assert_eq!(bitstream.read_uvlc(), 4);
        assert_eq!(bitstream.read_uvlc(), 5);
        assert_eq!(bitstream.read_uvlc(), 6);
        assert_eq!(bitstream.read_uvlc(), 7);
        assert_eq!(bitstream.read_uvlc(), 8);
        assert_eq!(bitstream.read_uvlc(), 9);
        assert_eq!(bitstream.read_uvlc(), 10);
        assert_eq!(bitstream.read_uvlc(), 11);
        assert_eq!(bitstream.read_uvlc(), 12);
        assert_eq!(bitstream.read_uvlc(), 13);
        assert_eq!(bitstream.read_uvlc(), 14);
    }

    #[test]
    fn test_bitstream_read_svlc() {
        let bitstream = Bitstream::new(vec![
            0b10100110, 0b01000010, 0b10011000, 0b11100010, 0b00000100, 0b10001010, 0b00010110,
            0b00110000, 0b01101000, 0b11100001, 0b11100000,
        ]);
        assert_eq!(bitstream.read_svlc(), 0);
        assert_eq!(bitstream.read_svlc(), 1);
        assert_eq!(bitstream.read_svlc(), -1);
        assert_eq!(bitstream.read_svlc(), 2);
        assert_eq!(bitstream.read_svlc(), -2);
        assert_eq!(bitstream.read_svlc(), 3);
        assert_eq!(bitstream.read_svlc(), -3);
        assert_eq!(bitstream.read_svlc(), 4);
        assert_eq!(bitstream.read_svlc(), -4);
        assert_eq!(bitstream.read_svlc(), 5);
        assert_eq!(bitstream.read_svlc(), -5);
        assert_eq!(bitstream.read_svlc(), 6);
        assert_eq!(bitstream.read_svlc(), -6);
        assert_eq!(bitstream.read_svlc(), 7);
        assert_eq!(bitstream.read_svlc(), -7);
    }
    #[test]
    fn test_copy_from() {
        let mut bitstream = Bitstream::new(vec![0b10101010, 0b11110000, 0b11001001, 0b00110011]);
        let bitstream2 = Bitstream::new(vec![0b11001001, 0b00110011, 0b11001001, 0b11111111]);
        bitstream.copy_from(&bitstream2, 1, 2);
        assert_eq!(
            bitstream.data,
            vec![0b00110011, 0b11001001, 0b11001001, 0b00110011]
        );
        bitstream.copy_from(&bitstream2, 3, 1);
        assert_eq!(
            bitstream.data,
            vec![0b00110011, 0b11001001, 0b11111111, 0b00110011]
        );
        bitstream.copy_from(&bitstream2, 0, 4);
        assert_eq!(
            bitstream.data,
            vec![
                0b00110011, 0b11001001, 0b11111111, 0b11001001, 0b00110011, 0b11001001, 0b11111111
            ]
        );
    }
}
