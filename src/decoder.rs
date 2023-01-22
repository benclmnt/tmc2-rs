use crate::{
    bitstream::{
        reader::{
            NalUnitType, PatchDataUnit, PatchModeITile, PatchModePTile, SeiPayloadType, TileType,
        },
        VideoBitstream, VideoType,
    },
    codec::{self, generate_point_cloud, EomParams, GroupOfFrames, Point3D, PointSet3},
    common::{
        context::{AtlasContext, AtlasFrameContext, TileContext},
        ColorFormat, VideoAttribute, VideoGeometry, VideoOccupancyMap,
    },
};

use super::common::context::Context;
use log::{debug, trace};
use num_enum::FromPrimitive;
use std::path::PathBuf;
use std::{cell::RefCell, fs::File};
use std::{io::Write, marker::PhantomData};
#[derive(Debug, Default)]
pub struct Params {
    pub start_frame: usize,
    pub compressed_stream_path: PathBuf,
    pub reconstructed_data_path: PathBuf,
    pub video_decoder_path: Option<PathBuf>,
    // (2Jan23): always true
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
    pub fn new(compressed_stream: PathBuf, video_decoder_path: Option<PathBuf>) -> Self {
        Self {
            compressed_stream_path: compressed_stream.clone(),
            reconstructed_data_path: PathBuf::from(&format!(
                "{:?}_dec_%04d.ply",
                compressed_stream.file_stem().unwrap_or_default()
            )),
            video_decoder_path,
            ..Default::default()
        }
    }

    pub fn with_start_frame(mut self, start_frame: usize) -> Self {
        self.start_frame = start_frame;
        self
    }

    pub fn with_reconstructed_data_path(mut self, reconstructed_data_path: PathBuf) -> Self {
        self.reconstructed_data_path = reconstructed_data_path;
        self
    }
}

pub struct Decoder {
    pub params: Params,
}

impl Decoder {
    pub fn new(params: Params) -> Self {
        Self { params }
    }

    pub fn decode(&self, context: &mut Context) -> GroupOfFrames {
        // TODO(12Dec22): this function can be parallelized it seems
        // if ( params_.nbThread_ > 0 ) { tbb::task_scheduler_init init( static_cast<int>( params_.nbThread_ ) ); }

        let mut atlas_context = Self::create_patch_frame(context);

        let sps = context.get_vps().expect("VPS not found");
        let ai = &sps.attribute_information;
        let oi = &sps.occupancy_information;
        let gi = &sps.geometry_information;
        let asps = context.get_atlas_sequence_parameter_set(0);
        let frame_count = atlas_context.frame_count();
        let ptl = &sps.profile_tier_level;
        let geometry_bitdepth = gi.geometry_2d_bitdepth_minus1 + 1;

        // TODO: maybe move this to context struct?
        let has_aux_data = asps.raw_patch_enabled_flag
            && asps.auxiliary_video_enabled_flag
            && sps.auxiliary_video_present_flag;
        assert!(!has_aux_data);

        // TODO: set_consitant_four_cc_code(context, 0);
        let occupancy_codec_id = CodecId::from(oi.occupancy_codec_id);
        let geometry_codec_id = CodecId::from(gi.geometry_codec_id);
        let _path_prefix = format!(
            "{:?}_dec_GOF{}",
            &self
                .params
                .compressed_stream_path
                .file_stem()
                .unwrap_or_default(),
            sps.v3c_parameter_set_id,
        );

        debug!(
            "CodecCodecId: profile_codec_group_idc={}, occupancy={:?}, geo={:?}",
            ptl.profile_codec_group_idc, occupancy_codec_id, geometry_codec_id
        );

        let mut video_decoder = LibavcodecDecoder {};

        let occ_bitstream = context
            .get_video_bitstream(VideoType::Occupancy)
            .expect("No occupancy bitstream");
        debug!(
            "\n*******Video Decoding: Occupancy ******** size: {}",
            occ_bitstream.len()
        );
        let occ_video: VideoOccupancyMap = video_decoder
            .decompress(
                occ_bitstream,
                VideoDecoderOptions {
                    codec_id: occupancy_codec_id,
                    bytestream_video_coder: true,
                    output_bitdepth: 8,
                    keep_intermediate_files: self.params.keep_intermediate_files,
                },
            )
            .unwrap();
        atlas_context.occ_frames = occ_video;
        assert_eq!(oi.occupancy_2d_bitdepth_minus1, 7);
        assert!(!oi.occupancy_msb_align_flag);
        // context.getVideoOccupancyMap().convertBitdepth( 8, oi.getOccupancy2DBitdepthMinus1() + 1,
        //                                           oi.getOccupancyMSBAlignFlag() );

        if sps.multiple_map_streams_present_flag {
            unimplemented!("multiple map streams not implemented");
        } else {
            trace!("Geometry\nMapIdx = 0, AuxiliaryVideoFlag = 0\n");
            let geo_bitstream = context
                .get_video_bitstream(VideoType::Geometry)
                .expect("No geometry bitstream");
            debug!(
                "\n*******Video Decoding: Geometry ******** size: {}",
                geo_bitstream.len()
            );
            let geo_video: VideoGeometry = video_decoder
                .decompress(
                    geo_bitstream,
                    VideoDecoderOptions {
                        codec_id: geometry_codec_id,
                        bytestream_video_coder: true,
                        output_bitdepth: geometry_bitdepth,
                        keep_intermediate_files: self.params.keep_intermediate_files,
                    },
                )
                .unwrap();

            atlas_context.geo_frames.push(geo_video);
            assert!(!gi.geometry_msb_align_flag);
            // TODO: convert bitdepth
            // context.getVideoGeometryMultiple(0).convertBitdepth( geometry_bitdepth, gi.getGeometry2DBitdepthMinus1() + 1, gi.geometry_msb_align_flag );
        }

        if has_aux_data {
            unimplemented!("Auxiliary geometry video not implemented")
        }

        // We only have 1 attribute actually
        assert_eq!(ai.attribute_count, 1);
        for i in 0..ai.attribute_count {
            let attribute_bitdepth = ai.attribute_2d_bitdepth_minus1[i as usize] + 1;
            // unused
            let attribute_type_id = ai.attribute_type_id[i as usize];
            let attribute_partition_dimension =
                ai.attribute_dimension_partitions_minus1[i as usize] + 1;
            let attribute_codec_id = CodecId::from(ai.attribute_codec_id[i as usize]);
            assert_eq!(attribute_partition_dimension, 1);

            debug!("CodecId attributeCodecId = {:?}", attribute_codec_id);
            for j in 0..attribute_partition_dimension {
                if sps.multiple_map_streams_present_flag {
                    unimplemented!("multiple map streams not implemented");
                } else {
                    trace!(
                        "Attribute\nAttrIdx = {}, AttrPartIdx = {}, AttributeTypeId = {}, MapIdx = 0, AuxiliaryVideoFlag = 0\n",
                        i,
                        j,
                        attribute_type_id,
                    );
                    let attr_bitstream = context
                        .get_video_bitstream(VideoType::Attribute)
                        .expect("No attribute bitstream");
                    debug!(
                        "\n*******Video Decoding: Attribute ******** size: {}",
                        attr_bitstream.len()
                    );
                    let attr_video: VideoAttribute = video_decoder
                        .decompress(
                            attr_bitstream,
                            VideoDecoderOptions {
                                codec_id: attribute_codec_id,
                                bytestream_video_coder: true,
                                output_bitdepth: attribute_bitdepth,
                                keep_intermediate_files: self.params.keep_intermediate_files,
                            },
                        )
                        .unwrap();

                    atlas_context.attr_frames.push(attr_video);

                    if has_aux_data {
                        unimplemented!("Auxiliary attribute video not implemented")
                    }
                }
            }
        }

        // Reconstruction
        let mut gof = GroupOfFrames {
            frames: Vec::with_capacity(frame_count),
        };
        // TODO: recreating the prediction list per attribute (either the attribtue is coded absolute, or follows the geometry)

        debug!("generate point cloud of {} frames", frame_count);
        // TODO: Looks like this is embarassingly parallel
        assert!(frame_count <= atlas_context.frame_count());
        for frame_idx in 0..frame_count {
            // all video have been decoded, start reconstruction processes
            if has_aux_data {
                unimplemented!("Auxiliary data not implemented")
            }

            let occupancy_precision = sps.frame_width as usize / atlas_context.occ_frames.width();
            let mut reconstruct = PointSet3::default();
            if ai.attribute_count > 0 {
                reconstruct.add_colors();
            }

            assert_eq!(
                atlas_context
                    .get_frame_context(frame_idx)
                    .num_tiles_in_atlas_frame,
                1,
                "we only support 1 tile per frame for now"
            );
            let frame_context = atlas_context.get_frame_context(frame_idx);
            let num_tiles_in_atlas_frame = frame_context.num_tiles_in_atlas_frame as usize;
            std::mem::drop(frame_context);
            for tile_idx in 0..num_tiles_in_atlas_frame {
                let atgl_idx = context
                    .atlas_hls
                    .get_atlas_tile_layer_index(frame_idx, tile_idx);
                assert_eq!(
                    atgl_idx, 0,
                    "looks like that's the case for decoding after reading the code..."
                );
                let gpc_params =
                    self.new_generate_point_cloud_params(context, atgl_idx, occupancy_precision);
                // pp_sei_params is mostly the same as gpc_params with some difference in handling attribute smoothing,
                // some of which we will ignore as I don't quite understand it yet
                //
                // < is postProcessing, > is generate point cloud
                // 143c143
                // <         for ( size_t i = 0; i < sei->getInstancesUpdated( k ); i++ ) {
                // ---
                // >         for ( size_t i = 0; i < sei->getInstancesUpdated( k ) + 1; i++ ) {
                // 149c149
                // <
                // ---
                // >             if ( sei->getMethodType( k, m ) == 1 ) { <--- in AttributeSmoothing if-clause
                // 182a187,188
                // >   params.useAuxSeperateVideo_           = asps.getAuxiliaryVideoEnabledFlag();
                // >
                let mut pp_sei_params =
                    self.new_generate_point_cloud_params(context, atgl_idx, occupancy_precision);
                pp_sei_params.use_aux_separate_video = false;

                trace!("processing frame {}, tile {}", frame_idx, tile_idx);
                if pp_sei_params.pbf.is_some() {
                    unimplemented!("PBF not implemented")
                }

                if num_tiles_in_atlas_frame > 1 {
                    unimplemented!("Multiple tiles not implemented")
                } else {
                    let mut frame_context = atlas_context.get_mut_frame_context(frame_idx);
                    let tile = frame_context.get_tile_mut(tile_idx);
                    let block_to_patch = codec::generate_block_to_patch_from_occupancy_map_video(
                        tile,
                        &atlas_context.occ_frames.frames[frame_idx],
                        1 << asps.log2_patch_packing_block_size,
                        occupancy_precision,
                    );
                    tile.block_to_patch = block_to_patch;
                }

                let (tile_construct, _partition) = {
                    let mut frame_context = atlas_context.get_mut_frame_context(frame_idx);
                    let tile = frame_context.get_tile_mut(tile_idx);
                    tile.frame_index = frame_idx;
                    tile.tile_index = tile_idx;
                    generate_point_cloud(
                        context,
                        &atlas_context,
                        tile,
                        &gpc_params,
                        true,
                        ai.attribute_count,
                    )
                    .unwrap()
                };

                if num_tiles_in_atlas_frame > 1 {
                    unimplemented!("Multiple tiles not implemented")
                }

                reconstruct.append_point_set(tile_construct);
            }

            // post processing
            // NOTE: atgl_index is always 0 when decoding. See how get_atlas_tile_layer_index is implemented. :)
            let pp_sei_params =
                self.new_generate_point_cloud_params(context, 0, occupancy_precision);
            trace!(
                "Post-Processing: postprocessSmoothing={} pbfEnableFlag={}",
                self.params.attr_transfer_filter_type,
                pp_sei_params.pbf.is_some()
            );

            if self.params.apply_geo_smoothing_type && pp_sei_params.geometry_smoothing.is_some() {
                unimplemented!("Geometry smoothing not implemented")
            }

            if ai.attribute_count > 0 {
                if self.params.apply_attr_smoothing_type && pp_sei_params.color_smoothing.is_some()
                {
                    unimplemented!("Color smoothing not implemented")
                }

                match atlas_context.attr_frames.get(0).unwrap().color_format() {
                    ColorFormat::_Rgb444 => reconstruct.copy_rgb16_to_rgb8(),
                    ColorFormat::Unknown => unreachable!(),
                    _ => reconstruct.convert_yuv16_to_rgb8(),
                }
            }

            // SKIP: some checksum stuff

            gof.frames.push(reconstruct);
        }

        gof
    }

    /// Create the context for all patches in all frames in atlas_context
    fn create_patch_frame(context: &mut Context) -> AtlasContext {
        let mut atlas_ctx = AtlasContext {
            frame_contexts: Vec::with_capacity(context.atlas_tile_layer_len()),
            ..Default::default()
        };
        let mut frame_count = 0;
        Self::set_tile_partition_size_afti(context);

        for i in 0..context.atlas_tile_layer_len() {
            let (afoc_msb, afoc_val) = context.derive_afoc_val(i);
            let atgl = context.get_mut_atlas_tile_layer(i);
            atgl.atlas_frame_order_count_msb = afoc_msb;
            atgl.atlas_frame_order_count_val = afoc_val;
            atgl.header.frame_index = afoc_val as u8;
            frame_count = std::cmp::max(frame_count, afoc_val + 1);
        }

        // context.atlas_contexts.resize(frame_count);
        // set_point_local_reconstruction(context);

        for atgl_idx in 0..context.atlas_tile_layer_len() {
            // (15Dec22) hmm this looks like the if clause always evaluates to true.
            // In what condition will the afoc_val be the same for 2 consecutive atlas tile layer?
            let mut afc = if atgl_idx == 0
                || context
                    .get_atlas_tile_layer(atgl_idx)
                    .atlas_frame_order_count_val
                    != context
                        .get_atlas_tile_layer(atgl_idx - 1)
                        .atlas_frame_order_count_val
            {
                let atgl = context.get_atlas_tile_layer(atgl_idx);
                let frame_index = atgl.header.frame_index;
                let afps_id = atgl.header.atlas_frame_parameter_set_id;
                Self::set_tile_size_and_location(context, frame_index as usize, afps_id as usize)
            } else {
                unreachable!("Looks like the if-clause will always evaluate to true");
                // AtlasFrameContext::default()
            };

            trace!("create_patch_frame_data_structure tile {}", atgl_idx);
            let sps = context.get_vps().unwrap();
            let atlu = context.get_atlas_tile_layer(atgl_idx);
            let ath = &atlu.header;
            let afps =
                context.get_atlas_frame_parameter_set(ath.atlas_frame_parameter_set_id as usize);
            let asps = context
                .get_atlas_sequence_parameter_set(afps.atlas_sequence_parameter_set_id as usize);
            let afti = &afps.atlas_frame_tile_information;
            let atgdu = &atlu.data_unit;
            // let geometry_bitdepth_2d = asps.geometry_2d_bitdepth_minus1 + 1;
            let geometry_bitdepth_3d = asps.geometry_3d_bitdepth_minus1 + 1;
            let frame_index = ath.frame_index;
            let tile_index = if afti.signalled_tile_id_flag {
                // afti.getTileId( ath.id )
                0
            } else {
                ath.id
            };
            assert_eq!(tile_index, 0);
            trace!(
                "create_patch_frame frame = {}, tiles = {}, atlas_index = {}, atgl_index = {}",
                frame_index,
                tile_index,
                0,
                atgl_idx
            );

            let patch_count = atgdu.patch_information_data.len();
            afc.title_frame_context = TileContext {
                frame_index: frame_index as usize,
                _atlas_frame_order_count_val: atlu.atlas_frame_order_count_val,
                _atlas_frame_order_count_msb: atlu.atlas_frame_order_count_msb,
                tile_index: tile_index as usize,
                _atl_index: atgl_idx,
                use_raw_points_separate_video: sps.auxiliary_video_present_flag
                    && asps.auxiliary_video_enabled_flag,
                _raw_patch_enabled_flag: asps.raw_patch_enabled_flag,
                _log2_patch_quantizer_size: ath.patch_size_info_quantizer,
                patches: Vec::with_capacity(patch_count),
                ..afc.title_frame_context
            };

            if frame_index > 0 && ath.tile_type != TileType::I {
                // (15Dec22) Observation from a few runs on the reference encoder.
                // IDEA: if there is only I tile, does it mean we can parallelize?
                unimplemented!("only I-tile in bitstream");
            }

            let tile_type = ath.tile_type;
            let min_level: usize = 1 << ath.pos_min_d_quantizer;
            // skipped: find the number of raw patches and eom patches

            // the for loop below is originally createPatchFrameDataStructure(context, atglIndex)
            // it inserts all patches to the frame context
            for patch_idx in 0..patch_count {
                let pid = &atgdu.patch_information_data[patch_idx];
                let patch_type =
                    PatchType::from_tile_type_and_patch_mode(tile_type, pid.patch_mode);

                match patch_type {
                    PatchType::Intra => {
                        let pdu = if let PatchDataUnit::Intra(pdu) = &pid.patch_data_unit {
                            pdu
                        } else {
                            unreachable!()
                        };
                        let packing_block_size = 1 << asps.log2_patch_packing_block_size;
                        let mut patch = Patch {
                            occupancy_resolution: packing_block_size,
                            uv0: pdu.pos_2d,
                            uv1: pdu.pos_3d_offset,
                            level_of_detail: if pdu.lod_enabled_flag {
                                unimplemented!("lod_enabled_flag = true not implemented");
                            } else {
                                (1, 1)
                            },
                            size_d: if pdu.pos_3d_range_d == 0 {
                                0
                            } else {
                                pdu.pos_3d_range_d * min_level - 1
                            },
                            size_uv0: if asps.patch_size_quantizer_present_flag {
                                (
                                    ((pdu.size_2d_minus1.0 + 1) as f64
                                        * (1 << ath.patch_size_info_quantizer.0) as f64
                                        / packing_block_size as f64)
                                        .ceil() as usize,
                                    ((pdu.size_2d_minus1.1 + 1) as f64
                                        * (1 << ath.patch_size_info_quantizer.1) as f64
                                        / packing_block_size as f64)
                                        .ceil() as usize,
                                )
                            } else {
                                (pdu.size_2d_minus1.0 + 1, pdu.size_2d_minus1.1 + 1)
                            },
                            _size_2d_in_pixel: if asps.patch_size_quantizer_present_flag {
                                (
                                    pdu.size_2d_minus1.0 * (1 << ath.patch_size_info_quantizer.0),
                                    pdu.size_2d_minus1.1 * (1 << ath.patch_size_info_quantizer.1),
                                )
                            } else {
                                (0, 0)
                            },
                            patch_orientation: pdu.orientation_index,
                            ..Default::default()
                        };
                        patch.set_view_id(pdu.projection_id);
                        if patch.projection_mode == 0 {
                            patch.d1 = pdu.pos_3d_offset_d * min_level;
                        } else {
                            let max_3d_coordinate = 1 << geometry_bitdepth_3d;
                            patch.d1 = max_3d_coordinate - (pdu.pos_3d_offset_d * min_level);
                        }
                        assert!(
                            patch.axes == (0, 2, 1)
                                || patch.axes == (1, 2, 0)
                                || patch.axes == (2, 0, 1)
                        );

                        trace!("patch(Intra) {}: UV0 {:?} UV1 {:?} D1={} S={:?} {}({}) P={} O={:?} A={:?} 45={} ProjId={} Axis={}", patch_idx, patch.uv0, patch.uv1, patch.d1, patch.size_uv1, patch.size_d, pdu.pos_3d_range_d, patch.projection_mode, patch.patch_orientation, patch.axes, asps.extended_projection_enabled_flag, pdu.projection_id, patch.axis_of_additional_plane);
                        // patch.alloc_one_layer_data (for PLR)
                        if asps.plr_enabled_flag {
                            unimplemented!("plr data not supported");
                        }
                        afc.title_frame_context.patches.push(patch);
                    }
                    PatchType::Inter => {
                        unimplemented!("inter patch not implemented");
                    }
                    PatchType::Merge => {
                        unimplemented!("merge patch not implemented");
                    }
                    PatchType::Skip => {
                        unreachable!("skip patch should not be in bitstream");
                    }
                    PatchType::_Raw => {
                        unimplemented!("raw patch not implemented");
                    }
                    PatchType::_Eom => {
                        unimplemented!("eom patch not implemented");
                    }
                    PatchType::_End => {
                        break;
                    }
                    _ => {
                        unreachable!("unknown frame / patch type");
                    }
                }
            }

            std::mem::drop(afps);
            atlas_ctx.frame_contexts.push(RefCell::new(afc));
        }

        // skipped: create hash sei for the last tile in frames.
        atlas_ctx
    }

    /// Sets partition's size (height and width) in AFTI
    pub(crate) fn set_tile_partition_size_afti(context: &mut Context) {
        let n = context.atlas_hls.atlas_frame_parameter_set.len();
        for i in 0..n {
            let mut afps = context.get_mut_atlas_frame_parameter_set(i);
            let asps = context
                .get_atlas_sequence_parameter_set(afps.atlas_sequence_parameter_set_id as usize);
            let frame_width = asps.frame_width;
            let frame_height = asps.frame_height;
            // NOTE (12Dec22): sets the partition_width and partition_height to asps.frame_width, asps.frame_height
            // ignore the rest since afti.single_tile_in_atlas_frame is always true
            let afti = &mut afps.atlas_frame_tile_information;
            if afti.single_tile_in_atlas_frame_flag {
                afti.set_partition_height(0, frame_height);
                afti.set_partition_width(0, frame_width);
            } else {
                unimplemented!("multiple tiles in atlas frame not implemented");
            }
        }
    }

    /// Set's tile size (height and width) and location (top_left)
    ///
    /// (16Dec22) Although it seems redundant, this function will be useful once we implement multi-tile in a frame
    pub(crate) fn set_tile_size_and_location(
        context: &mut Context,
        _frame_index: usize,
        afps_index: usize,
    ) -> AtlasFrameContext {
        let afps = context.get_atlas_frame_parameter_set(afps_index);

        let asps =
            context.get_atlas_sequence_parameter_set(afps.atlas_sequence_parameter_set_id as usize);
        let asps_auxiliary_video_enabled_flag = asps.auxiliary_video_enabled_flag;
        let frame_height = asps.frame_height;
        let frame_width = asps.frame_width;

        let afti_single_tile_in_atlas_frame_flag = afps
            .atlas_frame_tile_information
            .single_tile_in_atlas_frame_flag;

        std::mem::drop(afps);
        assert!(afti_single_tile_in_atlas_frame_flag);
        // a bunch of copying partition fields from afti to AtlasFrameContext is skipped.

        let afc = if afti_single_tile_in_atlas_frame_flag {
            AtlasFrameContext {
                _frame_width: frame_width,
                _frame_height: frame_height,
                num_tiles_in_atlas_frame: 1,
                title_frame_context: TileContext {
                    width: frame_width,
                    height: frame_height,
                    ..Default::default()
                },
                ..Default::default()
            }
        } else {
            unimplemented!("Only single tile in atlas frame is supported");
            // sets value for tiles' left_top, width, height,
            // atlas' frame_width, frame_height
        };

        if asps_auxiliary_video_enabled_flag {
            unimplemented!("Auxiliary video not supported");
            // sets value for atlas' aux_video_width, aux_video_tile_left_top
        }
        afc
    }

    /// Create parameters to generate point cloud
    fn new_generate_point_cloud_params(
        &self,
        context: &Context,
        atgl_index: usize,
        occupancy_precision: usize,
    ) -> codec::GeneratePointCloudParams {
        let sps = context.get_vps().expect("VPS not found");
        let gi = &sps.geometry_information;
        let asps = context.get_atlas_sequence_parameter_set(0);
        let _ptl = &sps.profile_tier_level;

        let params = codec::GeneratePointCloudParams {
            occupancy_resolution: 1 << asps.log2_patch_packing_block_size,
            occupancy_precision,
            enable_size_quantization: asps.patch_size_quantizer_present_flag,
            absolute_d1: sps.map_count_minus1 == 0 || sps.map_absolute_coding_enable_flag[1],
            multiple_streams: sps.multiple_map_streams_present_flag,
            _surface_thickness: asps.vpcc_extension.surface_thickness_minus1 + 1,
            remove_duplicate_points: self.params.point_local_reconstruction_type
                && asps.plr_enabled_flag,
            map_count_minus1: sps.map_count_minus1,
            single_map_pixel_interleaving: self.params.pixel_deinterleaving_type
                && asps.pixel_deinterleaving_flag,
            use_additional_points_patch: self.params.reconstruct_raw_type
                && asps.raw_patch_enabled_flag,
            use_aux_separate_video: asps.auxiliary_video_enabled_flag,
            enhanced_occupancy_map: if self.params.reconstruction_eom_type
                && asps.eom_patch_enabled_flag
            {
                Some(EomParams {
                    _fix_bitcount: asps.eom_fix_bit_count_minus1 + 1,
                })
            } else {
                None
            },
            _geometry_bitdepth_3d: gi.geometry_3d_coordinates_bitdepth_minus1 + 1,
            ..Default::default()
        };
        assert!(!params.multiple_streams, "multiple streams not supported");

        if self.params.apply_geo_smoothing_type
            && context.is_sei_present(
                NalUnitType::PrefixESEI,
                SeiPayloadType::GeometrySmoothing,
                atgl_index,
            )
        {
            unimplemented!("geometry smoothing not implemented")
        }

        if self.params.apply_occupancy_synthesis_type
            && context.is_sei_present(
                NalUnitType::PrefixESEI,
                SeiPayloadType::OccupancySynthesis,
                atgl_index,
            )
        {
            unimplemented!("occupancy synthesis not implemented")
        }

        if self.params.apply_attr_smoothing_type
            && context.is_sei_present(
                NalUnitType::PrefixESEI,
                SeiPayloadType::AttributeSmoothing,
                atgl_index,
            )
        {
            unimplemented!("attribute smoothing not implemented")
        }

        params
    }
}

enum PatchType {
    Intra,
    Inter,
    Merge,
    Skip,
    _Raw,
    _Eom,
    _End,
    Error,
}

impl PatchType {
    fn from_tile_type_and_patch_mode(tile_type: TileType, patch_mode: u8) -> Self {
        match tile_type {
            TileType::Skip => PatchType::Skip,
            TileType::P => match PatchModePTile::from(patch_mode) {
                PatchModePTile::Intra => PatchType::Intra,
                PatchModePTile::Inter => PatchType::Inter,
                PatchModePTile::Merge => PatchType::Merge,
                PatchModePTile::Skip => PatchType::Skip,
                _ => PatchType::Error,
            },
            TileType::I => match PatchModeITile::from(patch_mode) {
                PatchModeITile::Intra => PatchType::Intra,
                _ => PatchType::Error,
            },
        }
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, FromPrimitive)]
#[repr(u8)]
pub(crate) enum PatchOrientation {
    #[default]
    Default,
    Swap,
    Rot90,
    Rot180,
    Rot270,
    Mirror,
    MRot90,
    MRot180,
    MRot270,
}

/// Originally PCCPatch
#[derive(Default, Clone)]
pub(crate) struct Patch {
    // patch_index: usize,
    // original_index: usize,
    // frame_index: usize,
    // tile_index: usize,
    // index_in_frame: usize,
    /// u1: tangential shift
    /// v1: bitangential shift
    pub(crate) uv1: (usize, usize),
    size_uv1: (usize, usize),
    /// d1: depth shift
    pub(crate) d1: usize,
    /// size for depth
    size_d: usize,
    // /// size D pixel
    // size_d_pixel: usize,
    /// location in packed image (n * occupancy_resolution)
    pub(crate) uv0: (usize, usize),
    /// size of occupancy map (n * occupancy resolution)
    pub(crate) size_uv0: (usize, usize),
    _size_2d_in_pixel: (usize, usize),
    pub(crate) occupancy_resolution: usize,

    /// (LodScaleX, LodScaleYIdc)
    level_of_detail: (usize, usize),
    /// 0: related to min depth value, 1: related to the max value
    pub(crate) projection_mode: u8,
    /// x: normal axis, y: tangent axis, z: bitangent axis
    pub(crate) axes: (u8, u8, u8),
    pub(crate) axis_of_additional_plane: u8,
    // depth: (i16, i16),
    // /// occupancy map
    // occupancy: Vec<bool>,
    // /// view_id in 0..=5
    // view_id: u8,
    // /// index of matched patch from pre-frame patch
    // best_match_idx: i32,
    // ref_atlas_frame_idx: usize,
    // pred_type: usize,
    // /// Enhance delta depth
    // depth_eom: Vec<i16>,
    // /// for surface separation
    // depth_0pc_idx: Vec<i64>,
    /// patch orientation in canvas atlas
    pub(crate) patch_orientation: PatchOrientation,
    // point_local_reconstruction_lvl: u8,
    // point_local_reconstruction_mode_by_patch: u8,
    // point_local_reconstruction_mode_by_block: Vec<u8>,

    // cur_gpa_patch_data: GPAPatchData,
    // pre_gpa_patch_data: GPAPatchData,
    // is_global_patch: bool,
    // d0_count: usize,
    // eom_count: usize,
    // eom_and_d1_count: usize,
    // patch_type: u8,
    // is_roi_patch: bool,
    // roi_index: usize,
    // /// patch index
    // index_copy: usize,
    // bounding_box: Int16Box3D,
    // /// size of the depth map (width + 2 border)
    // border: i16,
    // depth_map_width: i16,
    // depth_map_height: i16,

    // /// list of neighboring patches' index
    // neighboring_patches: Vec<usize>,
    // depth_map: Vec<i16>,
    // occupancy_map: Vec<u8>,
    // /// 3D points created from borders of the patch
    // border_points: Vec<Point3D>,
}

impl Patch {
    /// Sets view id, axes, projection mode
    #[inline]
    fn set_view_id(&mut self, view_id: u8) {
        match view_id {
            0 => self.set_axis(0, 0, 2, 1, 0),
            1 => self.set_axis(0, 1, 2, 0, 0),
            2 => self.set_axis(0, 2, 0, 1, 0),

            3 => self.set_axis(0, 0, 2, 1, 1),
            4 => self.set_axis(0, 1, 2, 0, 1),
            5 => self.set_axis(0, 2, 0, 1, 1),

            6 => self.set_axis(1, 0, 2, 1, 0),
            7 => self.set_axis(1, 2, 0, 1, 0),
            8 => self.set_axis(1, 0, 2, 1, 1),
            9 => self.set_axis(1, 2, 0, 1, 1),

            10 => self.set_axis(2, 2, 0, 1, 0),
            11 => self.set_axis(2, 1, 2, 0, 0),
            12 => self.set_axis(2, 2, 0, 1, 1),
            13 => self.set_axis(2, 1, 2, 0, 1),

            14 => self.set_axis(3, 1, 2, 0, 0),
            15 => self.set_axis(3, 0, 2, 1, 0),
            16 => self.set_axis(3, 1, 2, 0, 1),
            17 => self.set_axis(3, 0, 2, 1, 1),
            _ => unreachable!(),
        }
    }

    #[inline]
    fn set_axis(&mut self, additional_plane: u8, normal: u8, tangent: u8, bitangent: u8, mode: u8) {
        self.axis_of_additional_plane = additional_plane;
        self.axes = (normal, tangent, bitangent);
        self.projection_mode = mode;
    }

    /// returns the index of the canvas block
    ///
    /// Originally PCCPatch::patchBlock2CanvasBlock. Ignore the tile parameters
    #[inline]
    pub(crate) fn patch_block_to_canvas_block(
        &self,
        u_blk: usize,
        v_blk: usize,
        canvas_stride_block: usize,
        canvas_height_block: usize,
    ) -> usize {
        let (x, y) = self.patch_to_canvas_helper(u_blk, v_blk, 1);
        assert!(x < canvas_stride_block && y < canvas_height_block);
        y * canvas_stride_block + x
    }

    #[inline]
    pub(crate) fn patch_to_canvas(
        &self,
        u: usize,
        v: usize,
        canvas_stride: usize,
        canvas_height: usize,
    ) -> (usize, usize) {
        let (x, y) = self.patch_to_canvas_helper(u, v, self.occupancy_resolution);
        assert!(x < canvas_stride && y < canvas_height);
        (x, y)
    }

    #[inline]
    fn patch_to_canvas_helper(&self, u: usize, v: usize, resolution: usize) -> (usize, usize) {
        let (u0, v0) = (self.uv0.0 * resolution, self.uv0.1 * resolution);
        let (size_u0, size_v0) = self.size_uv0;
        match self.patch_orientation {
            PatchOrientation::Default => (u + u0, v + v0),
            PatchOrientation::Rot90 => (size_v0 - 1 - v + u0, u + v0),
            PatchOrientation::Rot180 => (size_u0 - 1 - u + u0, size_v0 - 1 - v + v0),
            PatchOrientation::Rot270 => (v + u0, size_u0 - 1 - u + v0),
            PatchOrientation::Mirror => (size_u0 - 1 - u + u0, v + v0),
            PatchOrientation::MRot90 => (size_v0 - 1 - v + u0, size_u0 - 1 - u + v0),
            PatchOrientation::MRot180 => (u + u0, size_v0 - 1 - v + v0),
            PatchOrientation::MRot270 => (v + u0, u + v0),
            PatchOrientation::Swap => (v + u0, u + v0), // swap axis
        }
    }

    /// Coordinates might be truncated values.
    #[inline]
    pub(crate) fn generate_point(&self, u: usize, v: usize, depth: u16) -> Point3D {
        let mut point = Point3D::new(0, 0, 0);
        let (normal, tangent, bitangent) = self.axes;
        point[normal as usize] = self.generate_normal_coordinate(depth) as u16;
        point[tangent as usize] = (u * self.level_of_detail.0 + self.uv1.0) as u16;
        point[bitangent as usize] = (v * self.level_of_detail.1 + self.uv1.1) as u16;
        point
    }

    #[inline]
    pub(crate) fn generate_normal_coordinate(&self, depth: u16) -> usize {
        let depth = depth as usize;
        match self.projection_mode {
            0 => depth + self.d1,
            1 => std::cmp::max(self.d1, depth) - depth,
            _ => unreachable!(),
        }
    }
}

#[derive(Default, Debug, Clone, Copy)]
pub enum CodecId {
    H264,
    #[default]
    H265,
    H266,
}

impl CodecId {
    pub fn from(_codec_id: u8) -> CodecId {
        assert_eq!(_codec_id, 1);
        CodecId::H265
    }
}

#[derive(Default, Clone)]
pub(crate) struct Video<T>
where
    T: Copy + Default,
{
    frames: Vec<Image<T>>,
}

impl<T> Video<T>
where
    T: Copy + Default,
{
    #[inline]
    pub(crate) fn width(&self) -> usize {
        if self.frames.is_empty() {
            0
        } else {
            self.frames[0].width
        }
    }

    #[inline]
    pub(crate) fn height(&self) -> usize {
        if self.frames.is_empty() {
            0
        } else {
            self.frames[0].height
        }
    }

    #[inline]
    pub(crate) fn frame_count(&self) -> usize {
        self.frames.len()
    }

    #[inline]
    pub(crate) fn color_format(&self) -> ColorFormat {
        if self.frames.is_empty() {
            ColorFormat::Unknown
        } else {
            self.frames[0].format
        }
    }

    pub(crate) fn get(&self, index: usize) -> Option<&Image<T>> {
        self.frames.get(index)
    }
}

#[derive(Default, Clone, Debug)]
pub(crate) struct Image<T> {
    // can move this field to Video<T>?
    width: usize,
    height: usize,
    channels: [Vec<u8>; 3],
    format: ColorFormat,
    _phantom: PhantomData<T>,
}

impl<T> Image<T> {
    /// This is an internal function to get the index of the pixel in the image
    fn get_helper(&self, channel: usize, u: usize, v: usize) -> usize {
        assert!(channel < 3 && u < self.width as usize && v < self.height as usize);
        match (self.format, channel) {
            (ColorFormat::Yuv420, 0) => v * self.width + u,
            (ColorFormat::Yuv420, _) => (v / 2) * (self.width / 2) + (u / 2),
            _ => v * self.width + u,
        }
    }

    /// This is an internal function for debugging purposes. Please don't rely on it.
    #[allow(dead_code)]
    fn write(&self, filename: &str) {
        let mut file = File::create(filename).unwrap();
        let mut buf = Vec::new();
        // file.write_all(
        //     format!(
        //         "YUV4MPEG2 W{} H{} F30:1 Ip A0:0 C420mpeg2 XYSCSS=420MPEG2\n",
        //         self.width, self.height
        //     )
        //     .as_bytes(),
        // )
        // .unwrap();
        for channel in 0..3 {
            buf.extend_from_slice(&self.channels[channel]);
        }
        file.write_all(&buf).unwrap();
    }
}

impl Image<u8> {
    pub(crate) fn get(&self, channel: usize, u: usize, v: usize) -> u8 {
        let idx = self.get_helper(channel, u, v);
        self.channels[channel][idx]
    }
}

impl Image<u16> {
    /// We store image data (channels) internally as u8 slice. However, we want to be able to access as if it is stored as u16 slice.
    /// This function is used to achieve that.
    pub(crate) fn get(&self, channel: usize, u: usize, v: usize) -> u16 {
        let idx = self.get_helper(channel, u, v);
        // NOTE(21Jan23): we use from_ne_bytes here because libavcodec (the video decoder we use from ffmpeg) decodes to native endian
        u16::from_ne_bytes(
            self.channels[channel][2 * idx..2 * idx + 2]
                .try_into()
                .unwrap(),
        )
    }
}

#[derive(Default)]
struct VideoDecoderOptions {
    codec_id: CodecId,
    bytestream_video_coder: bool,
    output_bitdepth: u8,
    keep_intermediate_files: bool,
    // patch_color_subsampling: bool,
    // inverse_color_space_config: Option<PathBuf>,
    // color_space_conversion_path: Option<PathBuf>,
    // upsampling_filter: usize,
}

trait VideoDecoder {
    /// Does the heavylifting of calling the external video decoder (e.g. ffmpeg/HM/JM) and returns the decoded video.
    fn decode<T>(&mut self, data: &[u8], codec_id: CodecId) -> Result<Video<T>, ()>
    where
        T: Copy + Default;

    fn decompress<T>(
        &mut self,
        bitstream: &VideoBitstream,
        opts: VideoDecoderOptions,
    ) -> Result<Video<T>, ()>
    where
        T: Copy + Default,
    {
        // TODO: construct the filename for the decoded binstream
        let data = if opts.bytestream_video_coder {
            bitstream.sample_stream_to_bytestream(opts.codec_id, 4)
        } else {
            // TODO: optimize this expensive clone.
            bitstream.data.clone()
        };

        let video = self.decode(&data, opts.codec_id)?;
        // skipping some MD5 computation
        debug!(
            "Decoded video = {}x{} ({} frames) bitdepth={} color={:?}",
            video.width(),
            video.height(),
            video.frame_count(),
            opts.output_bitdepth,
            video.color_format()
        );

        // TODO(21Dec22): put human-distinguishable names for the files
        if opts.keep_intermediate_files {}

        // TODO(21Dec22): color conversion to yuv444

        // if opts.color_space_conversion_path.is_some() {
        //     unimplemented!()
        // } else {
        //     unimplemented!()
        // }

        // if opts.inverse_color_space_config.is_none() || video.is444() {
        // }
        Ok(video)
    }
}

pub struct LibavcodecDecoder {}

impl VideoDecoder for LibavcodecDecoder {
    fn decode<T>(&mut self, data: &[u8], codec_id: CodecId) -> Result<Video<T>, ()>
    where
        T: Copy + Default,
    {
        extern crate ffmpeg_next as ffmpeg;
        use ffmpeg::{codec, codec::context::Context, decoder, format, frame, Packet};
        use tempfile::NamedTempFile;

        let mut tmpfile = NamedTempFile::new().unwrap();
        println!("{:?}. bytestream length: {}", &tmpfile.path(), data.len());
        // TODO: use a buffer instead of writing to a tmpfile
        // Currently we write to disk because I can't figure out how to use ffmpeg API with a buffer
        // The tmpfile is not significant. It will be deleted when it goes out of scope
        tmpfile.write_all(data).unwrap();
        let mut ictx = format::input(&tmpfile.path()).unwrap();

        // transform to ffmpeg's codec id
        let codec = match codec_id {
            CodecId::H264 => codec::Id::H264,
            CodecId::H265 => codec::Id::HEVC,
            CodecId::H266 => codec::Id::VVC,
        };

        let mut decctx = Context::new();
        unsafe {
            let ptr = decctx.as_mut_ptr();
            (*ptr).pix_fmt = format::Pixel::YUV420P.into();
        }

        let mut decoder = decctx
            .decoder()
            .open_as(decoder::find(codec))
            .unwrap()
            .video()
            .unwrap();

        let mut video = Video::<T>::default();

        let mut frame_index = 0;
        let mut process_nal_unit = |packet: &Packet| {
            decoder.send_packet(packet).unwrap();
            let mut frame = frame::Video::empty();
            while decoder.receive_frame(&mut frame).is_ok() {
                let image = Image::<T> {
                    width: frame.width() as usize,
                    height: frame.height() as usize,
                    channels: [
                        frame.data(0).to_owned(),
                        frame.data(1).to_owned(),
                        frame.data(2).to_owned(),
                    ],
                    ..Default::default()
                };

                video.frames.push(image);
                frame_index += 1;
            }
        };

        for (_stream, packet) in ictx.packets() {
            // debug!("packet {}", _stream.index());
            process_nal_unit(&packet);
        }

        tmpfile.close().unwrap();
        Ok(video)
    }
}
