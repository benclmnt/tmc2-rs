# tmc2-rs

A port of [mpeg-pcc-tmc2](https://github.com/MPEGGroup/mpeg-pcc-tmc2). As of Jan 2023, it is based off [version 18](https://github.com/MPEGGroup/mpeg-pcc-tmc2/commit/30ce27aeb9d88d5dd2efedc6d45086396260444d).

## Features

- single-threaded (but at least 15x faster than the multi-threaded reference encoder)
- streaming library decoder API 
- (Jan 2023) Support rec0 profile_reconstruction_idc (this is a conformance level for the decoder)

## Building

1. Make sure you have ffmpeg 4.4 installed in your path. If you use nix, simply do `nix-shell`. Note that ffmpeg 5.x is NOT supported.
2. Run `cargo build --release --bin decoder`

## Generating bin files

To generate the bin files for testing,
1. Clone the [reference impl](https://github.com/MPEGGroup/mpeg-pcc-tmc2) and download a [dataset](http://plenodb.jpeg.org/pc/8ilabs/). 
2. Change `profileReconstructionIdc` at `cfg/common/ctc-common.cfg` to 0. As of Jan 2023, the decoder only supports rec0 profile_reconstruction_idc.

## Binary Usage

Simplest usage: `cargo run --bin decoder -- -i <COMPRESSED_STREAM_PATH> -o <RECONSTRUCTED_DATA_PATH>`

- `COMPRESSED_STREAM_PATH` is the path to your compressed bitstream input
- `RECONSTRUCTED_DATA_PATH` is the folder path for your decoded pointcloud.

See more options by running with `--help` flag.

### Env Variables

There are a few environment variables that can be set for better flexibility.
```
// show backtrace
export RUST_BACKTRACE=1 

// show log with granularity. By default only info-level (and above) logs are shown. 
// If the following value is exported, the decoder will show debug logs for all the modules in the library (except codec module), and show trace logs for the codec module.
export RUST_LOG=debug,tmc2rs::codec=trace 
```

## Library Usage

```rust
use std::path::PathBuf;
use tmc2rs::{Decoder, Params};

// Create the decoder, supplying the path to the compressed binstream
let mut decoder = Decoder::new(Params::new(PathBuf::from("path/to/compressed_stream")));

// Start the decoding process in a separate thread
decoder.start();

// Iterate through the decoded point cloud frames 
for frame in decoder.into_iter() {
   // do something with the frame
}
```
