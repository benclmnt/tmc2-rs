use clap::Parser;
use std::path::PathBuf;
use tmc2rs::bitstream;
use tmc2rs::bitstream::Bitstream;
use tmc2rs::common::context::Context;

/// An MPEG-VPCC-TMC2 conformant decoder
#[derive(Parser)]
struct Args {
    /// Configuration file name
    #[clap(short = 'c', long)]
    config: Option<String>,

    /// Output(encoder) / Input(decoder) compressed bitstream
    #[clap(long)]
    compressed_stream_path: PathBuf,

    /// Output decoded pointcloud. Multi-frame sequences maybe represented by %04i
    #[clap(long)]
    reconstructed_data_path: Option<PathBuf>,

    /// First frame number in sequence to encode/decode
    #[clap(long, default_value_t = 0)]
    start_frame: usize,

    /// Number of thread used for parallel processing
    #[clap(long, default_value_t = 4)]
    num_threads: u8,
    // missing: colorTransform, keepIntermediateFiles, colorSpaceConversionPath, inverseColorSpaceConversionConfig, videoDecoder{Occupancy,Geometry,Attribute}Path, patchColorSubsampling, shvcLayerIndex
}

impl Args {
    fn check(&self) -> bool {
        println!("Info: Using internal color space conversion");

        return true;
    }
}

fn decompress_video(args: Args) {
    let bitstream = Bitstream::from_file(&args.compressed_stream_path);
    let mut bitstream_stat = bitstream::Stat::new();
    // bitstream.computeMD5()
    // bitstreamStat.header = bitstream.size()
    let frame_number = args.start_frame;
    let (mut ssvu, header_size) =
        bitstream::reader::SampleStreamV3CUnit::from_bitstream(&bitstream);
    // bitstream_stat.incr_header(header_size);
    // DIFF: This is different (I think) from the reference implementation.
    let mut context = Context::new();
    while ssvu.get_v3c_unit_count() > 0 {
        // context.set_bitstream_stat(&bitstream_stat);
        ssvu.decode(&mut context);
        // TODO: context.check_profile()
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
