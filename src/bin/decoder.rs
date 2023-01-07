use clap::Parser;
use std::path::PathBuf;
use tmc2rs::bitstream;
use tmc2rs::bitstream::Bitstream;
use tmc2rs::common::context::Context;
use tmc2rs::decoder;

/// An MPEG-VPCC-TMC2 conformant decoder
#[derive(Parser)]
struct Args {
    /// Configuration file name
    #[clap(short = 'c', long)]
    config: Option<String>,

    /// Output(encoder) / Input(decoder) compressed bitstream
    #[clap(short = 'i', long)]
    compressed_stream_path: PathBuf,

    /// Output decoded pointcloud. Multi-frame sequences maybe represented by %04i
    #[clap(short = 'o', long)]
    reconstructed_data_path: Option<PathBuf>,

    /// First frame number in sequence to encode/decode
    #[clap(long, default_value_t = 0)]
    start_frame: usize,

    /// Number of thread used for parallel processing
    #[clap(long, default_value_t = 4)]
    num_threads: u8,

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
        println!("Info: Using internal color space conversion");

        true
    }
}

fn decompress_video(args: Args) {
    let bitstream = Bitstream::from_file(&args.compressed_stream_path);
    // let mut bitstream_stat = bitstream::Stat::new();
    // TODO[checks] bitstream.computeMD5()
    // TODO[stat] (9Dec22): Do everything related to bitstream_stat
    // bitstream_stat.header = bitstream.size()
    let frame_number = args.start_frame;
    let (mut ssvu, header_size) =
        bitstream::reader::SampleStreamV3CUnit::from_bitstream(&bitstream);
    // TODO[stat] bitstream_stat.incr_header(header_size);

    let decoder_params = decoder::Params::new(args.compressed_stream_path, args.video_decoder_path)
        .with_start_frame(args.start_frame);
    let decoder = decoder::Decoder::new(decoder_params);

    // IDEA (9Dec22): We can parallelize iterations of this loop, since the data is self-contained.
    // i.e. AD, OVD, GVD, AVD are independent only of the VPS that immediately precedes it.
    // In the reference implementation, after running `ssvu.decode(...)`, the decoder is run, which kinda implies that there is some potential for parallelism here.
    // Check how `context.active_vps` is updated.
    while ssvu.get_v3c_unit_count() > 0 {
        // DIFF: This is different (I think) from the reference implementation.
        let mut context = Context::default();
        // TODO[stat] context.set_bitstream_stat(&bitstream_stat);
        ssvu.decode(&mut context);
        // TODO[checks]: context.check_profile()

        // context.atlas_contexts[i].allocate_video_frames(&mut context);
        // context.atlas_index = atl_id as u8;
        let gof = decoder.decode(&mut context);
        // SKIP: a bunch of if caluses on metrics.
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

    decompress_video(args);
}
