use std::cell::{Ref, RefCell, RefMut};
use std::collections::hash_map::HashMap;
use std::ops::{Deref, DerefMut};
use std::rc::Rc;

use crate::bitstream::{
    reader,
    reader::{
        AtlasFrameParameterSetRbsp, AtlasSequenceParameterSetRbsp, AtlasTileHeader,
        AtlasTileLayerRbsp, V3CParameterSet, V3CUnitHeader, V3CUnitType,
    },
    VideoBitstream, VideoType,
};
use crate::common::{VideoAttribute, VideoGeometry, VideoOccupancyMap};
use crate::decoder::{Image, Patch, Video};
use cgmath::Vector3;

pub type Context = RawContext;

#[derive(Default)]
pub struct RawContext {
    // PCCContext.h
    // pub model_origin: Vector3<f64>,
    // pub model_scale: f64,
    // pub atlas_index: u8,
    // NOTE: we don't put atlas_context together with RawContext to avoid mutability issues
    // pub(crate) atlas_contexts: AtlasContext,

    // PCCHighLevelSyntax.h
    // pub(crate) video_bitstream: HashMap<VideoType, VideoBitstream>,
    v3c_unit_headers: HashMap<V3CUnitType, V3CUnitHeader>,
    vpcc_parameter_sets: Option<V3CParameterSet>,
    // pub active_vps: u8,

    // occupancy_precision: u8,
    // log2_patch_quantizer_size_x: u8,
    // log2_patch_quantizer_size_y: u8,
    // enable_patch_size_quantization: bool,
    // prefilter_lossy_OM: bool,
    // offset_lossy_OM: usize,
    // geometry_3d_coordinates_bitdepth: usize,
    // single_layer_mode: bool,
    // pub bitstream_stat: Stat,
    pub(crate) atlas_hls: AtlasHighLevelSyntax,
}

impl RawContext {
    pub(crate) fn get_v3c_unit_header(&self, v3c_unit_type: &V3CUnitType) -> Option<V3CUnitHeader> {
        self.v3c_unit_headers.get(v3c_unit_type).cloned()
    }

    pub(crate) fn set_v3c_unit_header(
        &mut self,
        v3c_unit_type: V3CUnitType,
        v3c_unit_header: V3CUnitHeader,
    ) {
        self.v3c_unit_headers.insert(v3c_unit_type, v3c_unit_header);
    }

    /// Add a new video bitstream to the context.
    pub(crate) fn add_video_bitstream(&mut self, video_bitstream: VideoBitstream) {
        self.atlas_hls.video_bitstreams.push(video_bitstream)
        // TODO[stat]: self.bitstream_stat.setVideoBinSize( video_bitstream.video_type, video_bitstream.size() );
    }

    /// Get a video bitstream from the context
    pub(crate) fn get_video_bitstream(&self, video_type: VideoType) -> Option<&VideoBitstream> {
        self.atlas_hls.get_video_bitstream(video_type)
    }

    /// Add a new V3CParameterSet to the context.
    #[inline]
    pub(crate) fn add_v3c_parameter_set(&mut self, v3c_parameter_set: V3CParameterSet) {
        assert!(self.vpcc_parameter_sets.is_none());
        self.vpcc_parameter_sets = Some(v3c_parameter_set);
    }

    /// Add a new ASPS to the context.
    #[inline]
    pub(crate) fn add_atlas_sequence_parameter_set(&mut self, asps: AtlasSequenceParameterSetRbsp) {
        self.atlas_hls.atlas_sequence_parameter_set.push(asps);
    }

    #[inline]
    pub(crate) fn get_atlas_sequence_parameter_set(
        &self,
        set_id: usize,
    ) -> &AtlasSequenceParameterSetRbsp {
        &self.atlas_hls.atlas_sequence_parameter_set[set_id]
    }

    #[inline]
    pub(crate) fn add_atlas_frame_parameter_set(&mut self, afps: AtlasFrameParameterSetRbsp) {
        self.atlas_hls
            .atlas_frame_parameter_set
            .push(RefCell::new(afps));
    }

    #[inline]
    pub(crate) fn get_atlas_frame_parameter_set(
        &self,
        set_id: usize,
    ) -> impl Deref<Target = AtlasFrameParameterSetRbsp> + '_ {
        Ref::map(
            self.atlas_hls.atlas_frame_parameter_set[set_id].borrow(),
            |s| s,
        )
    }

    #[inline]
    pub(crate) fn get_mut_atlas_frame_parameter_set(
        &self,
        set_id: usize,
    ) -> impl DerefMut<Target = AtlasFrameParameterSetRbsp> + '_ {
        RefMut::map(
            self.atlas_hls.atlas_frame_parameter_set[set_id].borrow_mut(),
            |s| s,
        )
    }

    /// Add a new Atlas Tile Layer to the context.
    #[inline]
    pub(crate) fn add_atlas_tile_layer(&mut self, atgl: AtlasTileLayerRbsp) {
        self.atlas_hls.atlas_tile_layer.push(atgl);
    }

    #[inline]
    pub(crate) fn atlas_tile_layer_len(&self) -> usize {
        self.atlas_hls.atlas_tile_layer.len()
    }

    #[inline]
    pub(crate) fn get_atlas_tile_layer(&self, layer_id: usize) -> &AtlasTileLayerRbsp {
        &self.atlas_hls.atlas_tile_layer[layer_id]
    }

    #[inline]
    pub(crate) fn get_mut_atlas_tile_layer(&mut self, layer_id: usize) -> &mut AtlasTileLayerRbsp {
        &mut self.atlas_hls.atlas_tile_layer[layer_id]
    }

    /// 8.4.3.1 Derives Atlas Frame Order Count (AFOC) value.
    ///
    /// Returns (afoc_msb, afoc_lsb). Note that afoc_val is equivalent to afoc_lsb
    pub(crate) fn derive_afoc_val(&self, atgl_index: usize) -> (u32, u32) {
        let atgh = &self.get_atlas_tile_layer(atgl_index).header;
        let afoc_lsb = atgh.atlas_frame_order_count_lsb;

        if atgl_index == 0 {
            return (0, afoc_lsb);
        }

        let afps = self.get_atlas_frame_parameter_set(atgh.atlas_frame_parameter_set_id as usize);
        let asps =
            self.get_atlas_sequence_parameter_set(afps.atlas_sequence_parameter_set_id as usize);
        let max_afoc_lsb = 1u32 << (asps.log2_max_atlas_frame_order_cnt_lsb_minus_4 + 4);

        assert!(atgl_index > 0);
        let prev_afoc_lsb = self
            .get_atlas_tile_layer(atgl_index - 1)
            .header
            .atlas_frame_order_count_lsb;
        let prev_afoc_msb = self
            .get_atlas_tile_layer(atgl_index - 1)
            .atlas_frame_order_count_msb;

        let afoc_msb = if afoc_lsb < prev_afoc_lsb && prev_afoc_lsb - afoc_lsb >= max_afoc_lsb / 2 {
            prev_afoc_msb + max_afoc_lsb
        } else if afoc_lsb > prev_afoc_lsb && afoc_lsb - prev_afoc_lsb > max_afoc_lsb / 2 {
            prev_afoc_msb - max_afoc_lsb
        } else {
            prev_afoc_msb
        };
        (afoc_msb, afoc_msb + afoc_lsb)
    }

    #[inline]
    pub(crate) fn allocate_atlas_hls(&mut self, size: usize) {
        assert!(
            size == 1,
            "V3C only have a single atlas. So there is no need to have a vector for atlas_hls."
        );
        // self.atlas_hls.resize_with(size, AtlasHighLevelSyntax::new);
    }

    #[inline]
    pub fn get_vps(&self) -> Option<&V3CParameterSet> {
        self.vpcc_parameter_sets.as_ref()
    }

    #[inline]
    pub(crate) fn get_num_ref_idx_active(&self, ath: &reader::AtlasTileHeader) -> usize {
        self.atlas_hls.get_num_ref_idx_active(ath)
    }

    pub(crate) fn is_sei_present(
        &self,
        nal_unit_type: reader::NalUnitType,
        payload_type: reader::SeiPayloadType,
        atgl_index: usize,
    ) -> bool {
        let atgl = self.get_atlas_tile_layer(atgl_index);
        if (*atgl.sei)
            .as_ref()
            .map(|sei| sei.is_sei_present(nal_unit_type, payload_type))
            .unwrap_or_default()
        {
            return true;
        }

        for i in atgl_index - 1..=0 {
            let atgl = self.get_atlas_tile_layer(i);
            if (*atgl.sei)
                .as_ref()
                .map(|sei| sei.is_sei_present(nal_unit_type, payload_type))
                .unwrap_or_default()
            {
                return true;
            }
        }
        false
    }
}

#[derive(Default)]
pub(crate) struct AtlasHighLevelSyntax {
    video_bitstreams: Vec<VideoBitstream>,
    pub(crate) atlas_sequence_parameter_set: Vec<AtlasSequenceParameterSetRbsp>,
    pub(crate) atlas_frame_parameter_set: Vec<RefCell<AtlasFrameParameterSetRbsp>>,
    // ref_atlas_frame_list: Vec<Vec<i32>>,
    // max_num_ref_atlas_frame: usize,
    // point_local_reconstruction_mode: Vec<PointLocalReconstructionMode>,
    pub(crate) atlas_tile_layer: Vec<AtlasTileLayerRbsp>,
}

impl AtlasHighLevelSyntax {
    pub fn get_num_ref_idx_active(&self, ath: &AtlasTileHeader) -> usize {
        use reader::TileType;

        let afps =
            &self.atlas_frame_parameter_set[ath.atlas_frame_parameter_set_id as usize].borrow();
        match ath.tile_type {
            TileType::I => 0,
            TileType::P | TileType::Skip => {
                if ath.num_ref_idx_active_override_flag {
                    ath.num_ref_idx_active_minus1 as usize + 1
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
                        afps.num_ref_idx_default_active_minus1 + 1,
                    ) as usize
                }
            }
        }
    }

    #[inline]
    fn get_video_bitstream(&self, video_type: VideoType) -> Option<&VideoBitstream> {
        self.video_bitstreams
            .iter()
            .rev()
            .find(|v| v.video_type == video_type)
    }

    /// Get the index of the atlas tile layer with the given frame and tile.
    #[inline]
    pub(crate) fn get_atlas_tile_layer_index(
        &self,
        frame_index: usize,
        tile_index: usize,
    ) -> usize {
        self.atlas_tile_layer
            .iter()
            .enumerate()
            .find(|(_, atgl)| {
                atgl.enc_frame_index == frame_index && atgl.enc_tile_index == tile_index
            })
            .map_or(0, |(i, _)| i)
    }
}

/// Context for a collection of frames
#[derive(Default, Clone)]
pub(crate) struct AtlasContext {
    // atlas_index: usize, // always 0 for V3C

    // encoding-only parameter
    // log2_max_atlas_frame_order_cnt_lsb: usize,
    pub(crate) frame_contexts: Vec<RefCell<AtlasFrameContext>>,
    // frames_in_afps: Vec<(usize, usize)>,
    pub(crate) occ_frames: VideoOccupancyMap,
    // occ_bitdepth: Vec<usize>,
    // occ_width: Vec<usize>,
    // occ_height: Vec<usize>,
    pub(crate) geo_frames: Vec<VideoGeometry>,
    // geo_bitdepth: Vec<usize>,
    // geo_width: Vec<usize>,
    // geo_height: Vec<usize>,
    // geo_aux_frames: VideoGeometry,
    pub(crate) attr_frames: Vec<VideoAttribute>,
    // attr_bitdepth: Vec<usize>,
    // attr_width: Vec<usize>,
    // attr_height: Vec<usize>,
    // attr_aux_frames: VideoAttribute,

    // GPA-related functions
    // sub_context: Vec<SubContext>,
    // union_patch: Vec<UnionPatch>,
}

impl AtlasContext {
    pub(crate) fn frame_count(&self) -> usize {
        self.frame_contexts.len()
    }

    pub(crate) fn get_frame_context(
        &self,
        frame_index: usize,
    ) -> impl Deref<Target = AtlasFrameContext> + '_ {
        self.frame_contexts[frame_index].borrow()
    }

    pub(crate) fn get_mut_frame_context(
        &self,
        frame_index: usize,
    ) -> impl DerefMut<Target = AtlasFrameContext> + '_ {
        self.frame_contexts[frame_index].borrow_mut()
    }
}

/// Context for an atlas frame.
/// It contains context for each patch that makes up the frame.
#[derive(Clone)]
pub(crate) struct AtlasFrameContext {
    // frame_index: usize,
    pub frame_width: u32,
    pub frame_height: u32,
    pub num_tiles_in_atlas_frame: u16, // (13Dec22) u16 because it is always 1. Originally usize
    /// (12Dec22) Always true
    pub single_partition_per_tile: bool,
    // uniform_partition_spacing: bool,
    // num_partition_cols: usize,
    // num_partition_rows: usize,
    // signalled_tile_id: bool,

    // /// These 4 partition fields are copied from AFTI. The value is set in PCCDecoder::setTileSizeAndLocation but never read.
    // partition_width: Vec<usize>,
    // partition_height: Vec<usize>,
    // partition_pos_x: Vec<usize>,
    // partition_pos_y: Vec<usize>,

    // aux_video_width: usize,
    // aux_tile_height: Vec<usize>,
    // aux_tile_left_top_y: Vec<usize>,

    // tile_id: Vec<usize>,
    // tile_contexts: Vec<TileContext>,
    // partition_to_tile_map: Vec<usize>,
    /// since `single_partition_per_tile` is true, we can change this to a non-vec field.
    ///
    /// Originally: titleFrameContext
    pub title_frame_context: TileContext,
}

impl Default for AtlasFrameContext {
    fn default() -> Self {
        Self {
            // frame_index: 0,
            frame_width: 0,
            frame_height: 0,
            num_tiles_in_atlas_frame: 0,
            single_partition_per_tile: true,
            title_frame_context: TileContext::default(),
        }
    }
}

impl AtlasFrameContext {
    #[inline]
    pub(crate) fn get_tile_mut(&mut self, tile_index: usize) -> &mut TileContext {
        assert_eq!(tile_index, 0, "we only support 1 tile per frame for now");
        assert_eq!(self.num_tiles_in_atlas_frame, 1);
        &mut self.title_frame_context
    }
}

/// Originally PCCFrameContext
#[derive(Default, Clone)]
pub(crate) struct TileContext {
    pub frame_index: usize,
    pub tile_index: usize,
    pub atl_index: usize,
    pub num_matched_patches: usize,
    pub(crate) width: u32,
    pub(crate) height: u32,
    pub left_top_in_frame: (usize, usize),
    // number_of_raw_points_patches: usize,
    // total_number_of_raw_points: usize,
    // total_number_of_eom_points: usize,
    pub total_number_of_regular_points: usize,
    pub global_patch_count: usize,
    pub geometry_3d_coordinates_bitdepth: usize,
    pub point_local_reconstruction_number: usize,
    pub use_raw_points_separate_video: bool,
    pub raw_patch_enabled_flag: bool,
    pub geometry_2d_bitdepth: usize,
    pub max_depth: usize,
    pub atlas_frame_order_count_val: u32,
    pub atlas_frame_order_count_msb: u32,
    pub atlas_frame_order_count_lsb: usize,
    // pub ref_afoc_list: Vec<u32>,
    pub num_ref_idx_active: usize,
    pub best_ref_list_index_in_asps: usize,
    pub referred_tile: bool,
    pub log2_patch_quantizer_size: (u8, u8),
    pub point_to_pixel: Vec<Vector3<usize>>,
    pub block_to_patch: Vec<usize>,
    pub occupancy_map: Vec<u8>,
    pub full_occupancy_map: Vec<u32>,
    pub patches: Vec<Patch>,
    // raw_points_patches: Vec<RawPointsPatch>,
    // num_of_raw_points: Vec<usize>,
    // raw_attributes: Vec<Color3B>,
    // eom_attributes: Vec<Color3B>,
    // src_point_cloud_by_patch: Vec<PointSet3>,
    // src_point_cloud_by_block: Vec<PointSet3>,
    // rec_point_cloud_by_block: Vec<PointSet3>,
    // point_to_pixel_by_block: Vec<Vec<Vector3<usize>>>,
    // pre_gpa_frame_size: GPAFrameSize,
    // cur_gpa_frame_size: GPAFrameSize,
    // ocp_gpa_info: FrameOcmInfo,
    // eom_patches: Vec<EomPatch>,
}
