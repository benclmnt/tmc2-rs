pub mod context;

pub(crate) enum CodecGroup {
    // AvcProgressiveHigh = 0,
    HevcMain10 = 1,
    Hevc444 = 2,
    // VvcMain10 = 3,
    // Mp4Ra = 127,
}
#[derive(Default, Debug, Clone, Copy)]
pub(crate) enum ColorFormat {
    Unknown,
    Rgb444,
    Yuv444,
    #[default]
    Yuv420,
}
