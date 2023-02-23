use clap::Parser;
use log::info;
use std::path::PathBuf;
use tmc2rs::writer::{Format, PlyWriter};
use tmc2rs::{Decoder, Params};

/// An MPEG-VPCC-TMC2 conformant decoder
#[derive(Parser)]
struct Args {
    // /// Configuration file name
    // #[clap(short = 'c', long)]
    // config: Option<String>,
    /// Path to the compressed bitstream input
    #[clap(short = 'i', long)]
    compressed_stream_path: PathBuf,

    /// Output path for decoded pointcloud.
    ///
    /// Specify either folder path or
    /// customized path for multi-frame sequences, with sequence number represented as %04d, e.g. `output/sequence_%04d.ply`
    #[clap(short = 'o', long)]
    reconstructed_data_path: PathBuf,

    /// First frame number in sequence to encode/decode
    #[clap(long, default_value_t = 0)]
    start_frame: usize,

    /// Number of thread used for parallel processing
    ///
    /// This option currently has no effect.
    #[clap(long, default_value_t = 4)]
    num_threads: u8,

    /// Keep intermediate files for inspection
    ///
    /// This option currently has no effect.
    #[clap(long, default_value_t = true)]
    keep_intermediate_files: bool,

    /// Path to the video decoder used to decompress occupancy, geometry, attribute maps
    #[clap(short = 'd', long)]
    video_decoder_path: Option<PathBuf>,
    // video_decoder_occupancy_path: PathBuf,
    // video_decoder_geometry_path: PathBuf,
    // video_decoder_attribute_path: PathBuf,

    // the following 3 attributes are always true
    // bytestream_video_coder_occupancy: bool,
    // bytestream_video_coder_geometry: bool,
    // bytestream_video_coder_attribute: bool,

    // /// The colour transform to be applied: (0): none, (1) RGB to YCbCr (Rec.709)
    // color_transform: ColorTransform,
    // /// Path to the HDRConvert. If unset, an internal color space onversion is used
    // color_space_conversion_path: PathBuf,
    // /// HDRConvert configuration file used for YUV420 to RGB444 conversion
    // inverse_color_space_conversion_config: PathBuf,

    // /// Enable per-patch color up-sampling
    // patch_color_subsampling: bool,

    // reconstruction options
    // pixel_deinterleaving_type: bool,
    // point_local_reconstruction_type: bool,
    // reconstruction_eom_type: bool,
    // duplicated_point_removal_type: bool,
    // reconstruct_raw_type: bool,
    // apply_geo_smooting_type: bool,
    // apply_attr_smoothing_type: bool,
    // attr_transfer_filter_type: bool,
    // apply_occupancy_synthesis_type: bool,
}

impl Args {
    fn check(&self) -> bool {
        // println!("Info: Using internal color space conversion");

        true
    }
}

fn main() {
    env_logger::init();

    println!("PccAppDecoder v18");
    let args: Args = Args::parse();
    if !args.check() {
        println!("Error: Input parameters are not correct");
        return;
    }

    let decoder_params = Params::new(args.compressed_stream_path);
    std::fs::create_dir_all(&args.reconstructed_data_path).unwrap();

    let mut decoder = Decoder::new(decoder_params);

    decoder.start();
    let path = args.reconstructed_data_path;

    for (i, frame) in decoder.into_iter().enumerate() {
        let frame_num = i + args.start_frame;
        let path = if path.is_dir() {
            path.join(format!("{frame_num:0>4}.ply"))
        } else {
            let parent = path.parent().unwrap();
            let filename = path.file_name().unwrap().to_str().unwrap();
            parent.join(filename.replace("%4d", format!("{frame_num:0>4}").as_ref()))
        };
        PlyWriter::new(frame, Format::Ascii).write(&path);
        info!("Frame {} written to {}", frame_num, path.display());
    }
}
