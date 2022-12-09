use std::collections::hash_map::HashMap;

use crate::bitstream::{
    reader,
    reader::{
        AtlasFrameParameterSetRbsp, AtlasSequenceParameterSetRbsp, AtlasTileLayerRbsp,
        V3CParameterSet, V3CUnitHeader, V3CUnitType,
    },
    AtlasHighLevelSyntax, Stat, VideoBitstream,
};
use cgmath::Vector3;

pub type Context = RawContext;
pub struct RawContext {
    // PCCContext.h
    pub model_origin: Vector3<f64>,
    pub model_scale: f64,
    // atlas_contexts: Vec<AtlasContext>,
    pub atlas_index: u8,

    // PCCHighLevelSyntax.h
    // video_bitstream: Vec<VideoBitstream>,
    v3c_unit_headers: HashMap<V3CUnitType, V3CUnitHeader>,
    vpcc_parameter_sets: Vec<V3CParameterSet>,
    pub active_vps: u8,
    // occupancy_precision: u8,
    // log2_patch_quantizer_size_x: u8,
    // log2_patch_quantizer_size_y: u8,
    // enable_patch_size_quantization: bool,
    // prefilter_lossy_OM: bool,
    // offset_lossy_OM: usize,
    // geometry_3d_coordinates_bitdepth: usize,
    // single_layer_mode: bool,
    // pub bitstream_stat: Stat,
    atlas_hls: Vec<AtlasHighLevelSyntax>,
}

impl RawContext {
    pub fn new() -> Self {
        Self {
            model_origin: Vector3::new(0.0, 0.0, 0.0),
            model_scale: 1.0,
            atlas_hls: Vec::new(),
            atlas_index: 0,
            active_vps: 0,
            v3c_unit_headers: HashMap::new(),
            vpcc_parameter_sets: Vec::new(),
        }
    }

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
        self.atlas_hls[self.atlas_index as usize]
            .video_bitstreams
            .push(video_bitstream)
        // TODO: self.bitstream_stat.setVideoBinSize( video_bitstream.video_type, video_bitstream.size() );
    }

    /// Add a new V3CParameterSet to the context.
    #[inline]
    pub(crate) fn add_v3c_parameter_set(&mut self, v3c_parameter_set: V3CParameterSet) {
        self.vpcc_parameter_sets.push(v3c_parameter_set);
    }

    /// Add a new ASPS to the context.
    #[inline]
    pub(crate) fn add_atlas_sequence_parameter_set(&mut self, asps: AtlasSequenceParameterSetRbsp) {
        self.atlas_hls[self.atlas_index as usize]
            .atlas_sequence_parameter_set
            .push(asps);
    }

    #[inline]
    pub(crate) fn get_atlas_sequence_parameter_set(
        &self,
        set_id: usize,
    ) -> &AtlasSequenceParameterSetRbsp {
        &self.atlas_hls[self.atlas_index as usize].atlas_sequence_parameter_set[set_id]
    }

    #[inline]
    pub(crate) fn add_atlas_frame_parameter_set(&mut self, afps: AtlasFrameParameterSetRbsp) {
        self.atlas_hls[self.atlas_index as usize]
            .atlas_frame_parameter_set
            .push(afps);
    }

    #[inline]
    pub(crate) fn get_atlas_frame_parameter_set(
        &self,
        set_id: usize,
    ) -> &AtlasFrameParameterSetRbsp {
        &self.atlas_hls[self.atlas_index as usize].atlas_frame_parameter_set[set_id]
    }

    /// Add a new Atlas Tile Layer to the context.
    #[inline]
    pub(crate) fn add_atlas_tile_layer(&mut self, atgl: AtlasTileLayerRbsp) {
        self.atlas_hls[self.atlas_index as usize]
            .atlas_tile_layer
            .push(atgl);
    }

    #[inline]
    pub(crate) fn allocate_atlas_hls(&mut self, size: usize) {
        self.atlas_hls.resize_with(size, AtlasHighLevelSyntax::new);
    }

    pub(crate) fn get_vps(&self) -> Option<&V3CParameterSet> {
        self.vpcc_parameter_sets.get(self.active_vps as usize)
    }

    pub(crate) fn get_num_ref_idx_active(&self, ath: &reader::AtlasTileHeader) -> usize {
        self.atlas_hls[self.atlas_index as usize].get_num_ref_idx_active(ath)
    }
}
