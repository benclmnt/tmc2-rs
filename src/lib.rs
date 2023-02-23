mod bitstream;
pub mod codec;
mod common;
mod decoder;
pub mod writer;

use bitstream::Bitstream;
use codec::PointSet3;
use common::context::Context;
use crossbeam_channel as chan;
use std::path::PathBuf;
use std::thread;

/// The library's decoder
pub struct Decoder {
    params: Params,
    // will be None once the decoder is started.
    tx: Option<chan::Sender<PointSet3>>,
    rx: chan::Receiver<PointSet3>,
}

/// Params to pass in to the library's decoder
#[derive(Debug, Default, Clone)]
pub struct Params {
    // NOTE: we don't need start_frame and reconstructed_data_path while decoding
    // pub start_frame: usize,
    // pub reconstructed_data_path: PathBuf,
    pub compressed_stream_path: PathBuf,
    pub video_decoder_path: Option<PathBuf>,
    // NOTE (2Jan23): always true
    // pub is_bytestream_video_coder: bool,
    pub keep_intermediate_files: bool,

    pub patch_color_subsampling: bool,
    pub color_space_conversion_path: Option<PathBuf>,
    pub inverse_color_space_conversion_config: Option<PathBuf>,

    // reconstruction options
    // NOTE (9Dec22): all set to default (false) for now since we are only supporting Rec0
    pixel_deinterleaving_type: bool,
    point_local_reconstruction_type: bool,
    reconstruction_eom_type: bool,
    _duplicated_point_removal_type: bool,
    reconstruct_raw_type: bool,
    apply_geo_smoothing_type: bool,
    apply_attr_smoothing_type: bool,
    attr_transfer_filter_type: bool,
    apply_occupancy_synthesis_type: bool,
}

impl Params {
    pub fn new(compressed_stream: PathBuf) -> Self {
        Self {
            compressed_stream_path: compressed_stream.clone(),
            ..Default::default()
        }
    }

    // pub fn with_start_frame(mut self, start_frame: usize) -> Self {
    //     self.start_frame = start_frame;
    //     self
    // }

    // pub fn with_video_decoder(mut self, video_decoder_path: PathBuf) -> Self {
    //     self.video_decoder_path = Some(video_decoder_path);
    //     self
    // }
}

impl Decoder {
    pub fn new(params: Params) -> Self {
        let (tx, rx) = chan::bounded(1);
        Self {
            params,
            tx: Some(tx),
            rx,
        }
    }

    /// Spawns a thread to decode.
    /// The decoded point cloud can be retrieved in order by repeatedly calling `recv_frame()` method until it returns None.
    ///
    /// Caller needs to ensure that this function is only called once per Decoder instance. Calling more than once will panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::path::PathBuf;
    /// use tmc2rs::{Decoder, Params};
    ///
    /// let mut decoder = Decoder::new(Params::new(PathBuf::from("path/to/compressed_stream"), None));
    /// decoder.start();
    /// for frame in decoder.into_iter() {
    ///    // do something with the frame
    /// }
    /// ```
    pub fn start(&mut self) {
        let bitstream = Bitstream::from_file(&self.params.compressed_stream_path);
        // let mut bitstream_stat = bitstream::Stat::new();
        // TODO[checks] bitstream.computeMD5()
        // TODO[stat] (9Dec22): Do everything related to bitstream_stat
        // bitstream_stat.header = bitstream.size()
        let (mut ssvu, _header_size) =
            bitstream::reader::SampleStreamV3CUnit::from_bitstream(&bitstream);
        // TODO[stat] bitstream_stat.incr_header(header_size);

        let decoder = decoder::Decoder::new(self.params.clone());
        let tx = self
            .tx
            .take()
            .expect("library decoder can only be started once");

        thread::spawn(move || {
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
                let _gof = decoder.decode(&mut context, tx.clone());
                // SKIP: a bunch of if clauses on metrics.
            }

            drop(tx);
        });
    }

    /// Blocks the current thread until the next decoded frame is received.
    ///
    /// Once this method returns None, it will not block anymore as there are no more frames left to be decoded.
    pub fn recv_frame(&self) -> Option<PointSet3> {
        self.rx.recv().ok()
    }
}

impl Iterator for Decoder {
    type Item = PointSet3;

    fn next(&mut self) -> Option<Self::Item> {
        self.recv_frame()
    }
}
