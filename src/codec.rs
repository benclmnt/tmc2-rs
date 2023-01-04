use std::path::PathBuf;

#[derive(Default, Clone, Debug)]
pub(crate) struct GeneratePointCloudParams {
    pub(crate) occupancy_resolution: usize,
    pub(crate) occupancy_precision: usize,
    pub(crate) enable_size_quantization: bool,
    pub(crate) neighbor_count_smoothing: usize,
    // radius2_smoothing: u64,
    // radius2_boundary_detection: u64,
    // pub(crate) raw_point_color_format: ColorFormat,
    // pub(crate) nb_thread: u8,
    pub(crate) multiple_streams: bool,
    pub(crate) absolute_d1: bool,
    pub(crate) surface_thickness: u8,
    // flagColorSmoothing is true if enhanced_occupancy_map is not None
    pub(crate) color_smoothing: Option<ColorSmoothingParams>,
    // flagGeometrySmoothing is true if enhanced_occupancy_map is not None
    pub(crate) geometry_smoothing: Option<GeometrySmoothingParams>,
    // enhancedOccupancyMapCode is true if enhanced_occupancy_map is not None
    pub(crate) enhanced_occupancy_map: Option<EomParams>,
    pub(crate) remove_duplicate_points: bool,
    pub(crate) map_count_minus1: u8,
    pub(crate) point_local_reconstruction: Option<PlrParams>,
    pub(crate) single_map_pixel_interleaving: bool,
    // pub(crate) path: PathBuf,
    pub(crate) use_additional_points_patch: bool,
    pub(crate) use_aux_separate_video: bool,
    pub(crate) geometry_bitdepth_3d: u8,
    // pub(crate) geometry_3d_coordinates_bitdepth: usize,
    /// Patch Block Filtering (pbf)
    pub(crate) pbf: Option<PbfParams>,
}

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct PbfParams {
    passes_count: i16,
    filter_size: i16,
    log2_threshold: i16,
    threshold_lossy_om: usize,
}

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct ColorSmoothingParams {
    cgrid_size: usize,
    threshold_color_smoothing: u64,
    threshold_color_difference: u64,
    threshold_color_variation: u64,
}

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct GeometrySmoothingParams {
    grid_size: usize,
    grid_smoothing: bool,
    threshold_smoothing: u64,
}

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct EomParams {
    pub(crate) fix_bitcount: u8,
}

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct PlrParams {
    pub(crate) number_of_modes: usize,
}
