pub mod context;

use crate::decoder::{Image, Video};

#[derive(Debug, Clone, Copy)]
pub(crate) enum CodecGroup {
    // AvcProgressiveHigh = 0,
    // HevcMain10 = 1,
    // Hevc444 = 2,
    // VvcMain10 = 3,
    // Mp4Ra = 127,
}

#[derive(Default, Debug, Clone, Copy)]
pub(crate) enum ColorFormat {
    Unknown,
    _Rgb444,
    // Yuv444,
    #[default]
    Yuv420,
}

pub(crate) type VideoOccupancyMap = Video<u8>;
pub(crate) type VideoGeometry = Video<u16>;
pub(crate) type VideoAttribute = Video<u16>;
pub(crate) type ImageOccupancyMap = Image<u8>;

pub(crate) const INTERMEDIATE_LAYER_INDEX: usize = 100;
