use log::debug;
use std::cell::RefCell;
use std::path::Path;

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
        return self.position.borrow().bytes;
    }

    /// Originally: Bitstream::capacity
    pub fn size(&self) -> u64 {
        return self.data.len() as u64;
    }

    pub fn more_data(&self) -> bool {
        return self.position_in_bytes() < self.data.len();
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
        return val;
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
        return (1 << leading_zeros) - 1 + self.read(leading_zeros);
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

pub(crate) struct AtlasHighLevelSyntax {
    pub video_bitstreams: Vec<VideoBitstream>,
    pub atlas_sequence_parameter_set: Vec<reader::AtlasSequenceParameterSetRbsp>,
    pub atlas_frame_parameter_set: Vec<reader::AtlasFrameParameterSetRbsp>,
    // ref_atlas_frame_list: Vec<Vec<i32>>,
    // max_num_ref_atlas_frame: usize,
    // point_local_reconstruction_mode: Vec<PointLocalReconstructionMode>,
    pub atlas_tile_layer: Vec<reader::AtlasTileLayerRbsp>,
}

impl AtlasHighLevelSyntax {
    pub fn new() -> Self {
        Self {
            video_bitstreams: Vec::new(),
            atlas_sequence_parameter_set: Vec::new(),
            atlas_frame_parameter_set: Vec::new(),
            atlas_tile_layer: Vec::new(),
        }
    }

    pub fn get_num_ref_idx_active(&self, ath: &reader::AtlasTileHeader) -> usize {
        use reader::TileType;

        let afps = &self.atlas_frame_parameter_set[ath.atlas_frame_parameter_set_id as usize];
        let num_ref_idx_active = match ath.tile_type {
            TileType::I => 0,
            TileType::P | TileType::Skip => {
                if ath.num_ref_idx_active_override_flag {
                    ath.num_ref_idx_active_minus_1 + 1
                } else {
                    let asps = &self.atlas_sequence_parameter_set
                        [afps.atlas_sequence_parameter_set_id as usize];
                    let ref_list = if ath.ref_atlas_frame_list_sps_flag {
                        &asps.ref_list_struct[ath.ref_atlas_frame_list_idx as usize]
                    } else {
                        &ath.ref_list_struct
                    };
                    std::cmp::min(
                        ref_list.num_ref_entries,
                        afps.num_ref_idx_default_active_minus_1 + 1,
                    )
                }
            }
        } as usize;
        num_ref_idx_active
    }
}

pub(crate) struct VideoBitstream {
    data: Vec<u8>,
    video_type: VideoType,
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
}

// TODO: set Attribute enum values
#[derive(Debug, FromPrimitive, IntoPrimitive)]
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
