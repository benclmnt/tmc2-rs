use std::{ops::Add, path::Path};

use crate::{
    common::{
        context::{AtlasContext, Context, TileContext},
        ImageOccupancyMap, VideoAttribute, VideoGeometry, INTERMEDIATE_LAYER_INDEX,
    },
    writer::{Format, PlyWriter},
};
use cgmath::Vector3;
use log::trace;

pub(crate) type Point3D = Vector3<u16>;
type Color3B = Vector3<u8>;

type Color16bit = Vector3<u16>;
// type Normal3D = Vector3<usize>;
// type Matrix3D = Matrix3<usize>;

#[derive(Debug, Default)]
pub(crate) struct PointSet3 {
    // NOTE: IF YOU UPDATE THIS STRUCT, dont forget to update resize and append point set.
    pub(crate) positions: Vec<Point3D>,
    pub(crate) colors: Vec<Color3B>,
    colors16bit: Vec<Color16bit>,
    // boundary_point_types: Vec<u16>,
    point_patch_indexes: Vec<(usize, usize)>,
    // parent_point_index: Vec<usize>,
    /// Only if PCC_SAVE_POINT_TYPE is true
    // types: Vec<u8>,
    // reflectances: Vec<u16>,
    // normals: Vec<Normal3D>,
    // with_normals: bool,
    pub(crate) with_colors: bool,
    // with_reflectances: bool,
}

impl PointSet3 {
    #[inline]
    pub(crate) fn point_count(&self) -> usize {
        self.positions.len()
    }

    /// add point to PointSet, and allocate the rest of the structure and returns its index
    pub(crate) fn add_point(&mut self, position: Point3D) -> usize {
        self.positions.push(position);
        if self.with_colors {
            self.colors.push(Color3B::new(127, 127, 127));
            self.colors16bit.push(Color16bit::new(0, 0, 0));
        }
        self.point_patch_indexes.push((0, 0));
        self.positions.len() - 1
    }

    /// set PointSet to use color
    pub(crate) fn add_colors(&mut self) {
        self.with_colors = true;
        self.reserve(self.point_count());
    }

    pub(crate) fn append_point_set(&mut self, pointset: PointSet3) -> usize {
        self.positions.extend(pointset.positions.iter());
        self.colors.extend(pointset.colors.iter());
        self.colors16bit.extend(pointset.colors16bit.iter());
        self.point_patch_indexes
            .extend(pointset.point_patch_indexes.iter());

        // SKIP: self.resize(self.point_count()). for what?
        self.point_count()
    }

    pub(crate) fn reserve(&mut self, size: usize) {
        self.positions.reserve(size);
        if self.with_colors {
            self.colors.reserve(size);
            self.colors16bit.reserve(size);
        }
        self.point_patch_indexes.reserve(size);
    }

    #[inline]
    /// add color to PointSet
    pub(crate) fn _set_color(&mut self, index: usize, color: Color3B) {
        assert!(self.with_colors && index < self.colors.len());
        self.colors[index] = color;
    }

    pub(crate) fn convert_yuv16_to_rgb8(&mut self) {
        assert!(self.with_colors);
        assert!(self.colors16bit.len() == self.point_count());
        for (i, color16bit) in self.colors16bit.iter().enumerate() {
            self.colors[i] = convert_yuv16_to_rgb8(color16bit);
        }
    }

    pub(crate) fn copy_rgb16_to_rgb8(&mut self) {
        assert!(self.with_colors);
        assert!(self.colors16bit.len() == self.point_count());
        for (i, color16bit) in self.colors16bit.iter().enumerate() {
            self.colors[i] = Vector3 {
                x: color16bit.x as u8,
                y: color16bit.y as u8,
                z: color16bit.z as u8,
            };
        }
    }
}

#[derive(Debug, Default)]
pub struct GroupOfFrames {
    pub(crate) frames: Vec<PointSet3>,
}

impl GroupOfFrames {
    pub fn load() -> bool {
        // TODO
        true
    }

    pub fn write(self, path: &Path) {
        for (i, frame) in self.frames.into_iter().enumerate() {
            let path = if path.is_dir() {
                path.join(format!("{:0>4}.ply", i))
            } else {
                let parent = path.parent().unwrap();
                let filename = path.file_name().unwrap().to_str().unwrap();
                parent.join(filename.replace("%4d", format!("{:0>4}", i).as_ref()))
            };
            PlyWriter::new(frame, Format::Ascii).write(&path);
        }
    }
}

#[derive(Default, Clone, Debug)]
pub(crate) struct GeneratePointCloudParams {
    pub(crate) occupancy_resolution: usize,
    pub(crate) occupancy_precision: usize,
    pub(crate) enable_size_quantization: bool,
    // pub(crate) neighbor_count_smoothing: usize,
    // radius2_smoothing: u64,
    // radius2_boundary_detection: u64,
    // pub(crate) raw_point_color_format: ColorFormat,
    // pub(crate) nb_thread: u8,
    pub(crate) multiple_streams: bool,
    pub(crate) absolute_d1: bool,
    pub(crate) _surface_thickness: u8,
    /// flagColorSmoothing is true if enhanced_occupancy_map is not None
    pub(crate) color_smoothing: Option<ColorSmoothingParams>,
    /// flagGeometrySmoothing is true if enhanced_occupancy_map is not None
    pub(crate) geometry_smoothing: Option<GeometrySmoothingParams>,
    /// enhancedOccupancyMapCode is true if enhanced_occupancy_map is not None
    pub(crate) enhanced_occupancy_map: Option<EomParams>,
    pub(crate) remove_duplicate_points: bool,
    pub(crate) map_count_minus1: u8,
    pub(crate) point_local_reconstruction: Option<PlrParams>,
    pub(crate) single_map_pixel_interleaving: bool,
    // pub(crate) path: PathBuf,
    pub(crate) use_additional_points_patch: bool,
    pub(crate) use_aux_separate_video: bool,
    pub(crate) _geometry_bitdepth_3d: u8,
    // pub(crate) geometry_3d_coordinates_bitdepth: usize,
    /// Patch Block Filtering (pbf)
    pub(crate) pbf: Option<PbfParams>,
}

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct PbfParams {
    _passes_count: i16,
    _filter_size: i16,
    _log2_threshold: i16,
    _threshold_lossy_om: usize,
}

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct ColorSmoothingParams {
    _cgrid_size: usize,
    _threshold_color_smoothing: u64,
    _threshold_color_difference: u64,
    _threshold_color_variation: u64,
}

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct GeometrySmoothingParams {
    _grid_size: usize,
    _grid_smoothing: bool,
    _threshold_smoothing: u64,
}

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct EomParams {
    pub(crate) _fix_bitcount: u8,
}

#[derive(Default, Debug, Clone, Copy)]
pub(crate) struct PlrParams {
    pub(crate) _number_of_modes: usize,
}

pub(crate) fn generate_block_to_patch_from_occupancy_map_video(
    tile: &TileContext,
    occupancy_map_image: &ImageOccupancyMap,
    occupancy_resolution: usize,
    occupancy_precision: usize,
) -> Vec<usize> {
    let patches = &tile.patches;
    let block_to_patch_width = tile.width / occupancy_resolution as u32;
    let block_to_patch_height = tile.height / occupancy_resolution as u32;
    let block_count = block_to_patch_width * block_to_patch_height;

    let mut block_to_patch: Vec<usize> = vec![0; block_count as usize];
    for (patch_index, patch) in patches.iter().enumerate() {
        for v0 in 0..patch.size_uv0.1 {
            for u0 in 0..patch.size_uv0.0 {
                let block_index = patch.patch_block_to_canvas_block(
                    u0,
                    v0,
                    block_to_patch_width as usize,
                    block_to_patch_height as usize,
                );
                let mut non_zero_pixel = 0usize;
                for v1 in 0..patch.occupancy_resolution {
                    let v = v0 * patch.occupancy_resolution + v1;
                    for u1 in 0..patch.occupancy_resolution {
                        let u = u0 * patch.occupancy_resolution + u1;
                        let (mut x, mut y) =
                            patch.patch_to_canvas(u, v, tile.width as usize, tile.height as usize);
                        x += tile.left_top_in_frame.0;
                        y += tile.left_top_in_frame.1;
                        non_zero_pixel += occupancy_map_image.get(
                            0,
                            x / occupancy_precision,
                            y / occupancy_precision,
                        ) as usize;
                    }
                }
                if non_zero_pixel > 0 {
                    block_to_patch[block_index] = patch_index + 1;
                }
            }
        }
    }

    block_to_patch
}

/// Returns a PointSet3 of reconstructed frame and a vector of partitions
///
/// DIFF: calls color_point_cloud here, instead asking the caller to do it.
/// So the resulting point cloud already contains attributes.
pub(crate) fn generate_point_cloud(
    context: &Context,
    atlas: &AtlasContext,
    tile: &mut TileContext,
    params: &GeneratePointCloudParams,
    _is_decoder: bool,
    attribute_count: u8,
) -> Option<(PointSet3, Vec<usize>)> {
    let frame_index = tile.frame_index;
    let tile_index = tile.tile_index;
    trace!("generate point cloud F = {} start", frame_index);
    assert!(!atlas.geo_frames.is_empty());
    let geo_video = &atlas.geo_frames[0];
    let block_to_patch_width = tile.width as usize / params.occupancy_resolution;
    let block_to_patch_height = tile.height as usize / params.occupancy_resolution;
    let map_count = params.map_count_minus1 as usize + 1;

    let mut reconstruct = PointSet3::default();
    if attribute_count > 0 {
        reconstruct.add_colors();
    }

    // most likely an overallocation but probably doesn't matter
    reconstruct.reserve(2 * tile.width as usize * tile.height as usize * tile.patches.len());
    let mut partition = vec![];

    trace!(
        "generate point cloud pbfEnableFlag = {} ",
        params.pbf.is_some()
    );

    if params.pbf.is_some() {
        unimplemented!("pbfEnableFlag is not implemented yet");
    } else {
        // point cloud occupancy map upscaling from video using nearest neighbour
        let height = tile.height as usize;
        let width = tile.width as usize;
        tile.occupancy_map = vec![0; height * width];
        for v in 0..height {
            for u in 0..width {
                tile.occupancy_map[v * width + u] = atlas.occ_frames.get(frame_index).unwrap().get(
                    0,
                    (tile.left_top_in_frame.0 + u) / params.occupancy_precision,
                    (tile.left_top_in_frame.1 + v) / params.occupancy_precision,
                );
            }
        }
    }

    if params.enable_size_quantization {
        unimplemented!("enableSizeQuantization is not implemented yet");
    }
    // SKIP: because we dont reuse struct like in the reference impl
    // point_to_pixel.resize(0)
    // reconstruct.clear()

    trace!("Frame {} in generatePointCloud, params.useAdditionalPointsPatch = {}, params.enhancedOccupancyMapCode = {}", frame_index, params.use_additional_points_patch, params.enhanced_occupancy_map.is_some());
    trace!("mapCount = {}", map_count);

    let video_frame_index;
    if params.multiple_streams {
        unimplemented!("multipleStreams is not implemented yet");
    } else {
        video_frame_index = tile.frame_index * map_count;
        if geo_video.frame_count() < video_frame_index + map_count {
            return None;
        }
    };

    trace!(
        "videoFrameIndex(shift):frameIndex*mapCount = {}",
        video_frame_index
    );
    let _frame0 = if params.multiple_streams {
        unimplemented!("multipleStreams is not implemented yet");
    } else {
        geo_video.get(video_frame_index).unwrap()
    };

    let _bp_flag = if params.pbf.is_none() {
        vec![0, tile.width as usize * tile.height as usize]
    } else {
        vec![]
    };
    // let mut eom_points_per_patch = vec![];
    // eom_points_per_patch.resize
    let patch_precedence_order_flag = context
        .get_atlas_sequence_parameter_set(0)
        .patch_precedence_order_flag;
    assert!(
        !patch_precedence_order_flag,
        "support for patchPrecedenceOrderFlag is not implemented yet"
    );
    // if decoder && patchPrecedenceOrderFlag, reverse the patch iterator

    // DIFF: we don't set this back to tile.point_to_pixel
    let mut point_to_pixel = vec![];

    for (patch_index, patch) in tile.patches.iter().enumerate() {
        trace!(
            "P{}/{}: 2D={:?}*{:?} 3D({},{},{})*({},{}) A={:?} Or={:?} P={} => {} AxisOfAdditionalPlane={}",
            patch_index,
            tile.patches.len(),
            patch.uv0,
            patch.size_uv0,
            patch.uv1.0,
            patch.uv1.1,
            patch.d1,
            patch.size_uv0.0 * patch.occupancy_resolution,
            patch.size_uv0.1 * patch.occupancy_resolution,
            patch.axes,
            patch.patch_orientation,
            patch.projection_mode,
            reconstruct.point_count(),
            patch.axis_of_additional_plane,
        );

        for v0 in 0..patch.size_uv0.1 {
            for u0 in 0..patch.size_uv0.0 {
                let block_index = patch.patch_block_to_canvas_block(
                    u0,
                    v0,
                    block_to_patch_width,
                    block_to_patch_height,
                );
                if tile.block_to_patch[block_index] != patch_index + 1 {
                    continue;
                }
                for v1 in 0..patch.occupancy_resolution {
                    let v = v0 * patch.occupancy_resolution + v1;
                    for u1 in 0..patch.occupancy_resolution {
                        let u = u0 * patch.occupancy_resolution + u1;
                        let (x, y) =
                            patch.patch_to_canvas(u, v, tile.width as usize, tile.height as usize);
                        let x_in_video_frame = x + tile.left_top_in_frame.0;
                        let y_in_video_frame = y + tile.left_top_in_frame.1;
                        if params.pbf.is_some() {
                            unimplemented!("pbfEnableFlag is not implemented yet");
                        } else {
                            let occupancy = tile.occupancy_map[y * tile.width as usize + x];
                            if occupancy == 0 {
                                continue;
                            }
                        }

                        if params.enhanced_occupancy_map.is_some() {
                            unimplemented!("enhancedOccupancyMapCode is not implemented yet");
                        } else {
                            let created_points = if params.point_local_reconstruction.is_some() {
                                unimplemented!("pointLocalReconstruction is not implemented yet");
                            } else {
                                generate_points(
                                    params,
                                    tile,
                                    geo_video,
                                    video_frame_index,
                                    patch_index,
                                    (u, v),
                                    (x_in_video_frame, y_in_video_frame),
                                )
                            };
                            assert_eq!(
                                created_points.len(),
                                2,
                                "looks like it from code comprehension"
                            );

                            for i in 0..created_points.len() {
                                if params.remove_duplicate_points
                                    && i == 0
                                    && created_points[i] == created_points[0]
                                {
                                    continue;
                                }
                                let _point_index = match patch.axis_of_additional_plane {
                                    0 => {
                                        let point_index = reconstruct.add_point(created_points[i]);
                                        reconstruct.point_patch_indexes[point_index] =
                                            (tile_index, patch_index);
                                        point_index
                                    }
                                    _ => unimplemented!(
                                        "axisOfAdditionalPlane is not implemented yet"
                                    ),
                                };

                                // NOTE(6Jan23): the reference implementation put color with random rgb here, but it is bogus
                                // the color will be updated by the caller of this fn afterwards, so there is no need to put color here
                                // reconstruct.set_color(point_index, color);

                                if params.pbf.is_some() {
                                    unimplemented!("pbf not implemented")
                                }

                                // SKIP: #cfg[PCC_SAVE_POINT_TYPE]

                                partition.push(patch_index);

                                if params.single_map_pixel_interleaving {
                                    unimplemented!(
                                        "singleMapPixelInterleaving is not implemented yet"
                                    );
                                } else if params.point_local_reconstruction.is_some() {
                                    unimplemented!(
                                        "pointLocalReconstruction is not implemented yet"
                                    );
                                } else {
                                    point_to_pixel.push(Vector3 {
                                        x,
                                        y,
                                        z: if i < 2 {
                                            i
                                        } else {
                                            // should be unreachable if created_points.len() is always 2.
                                            INTERMEDIATE_LAYER_INDEX + 1
                                        },
                                    })
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    tile.total_number_of_regular_points = reconstruct.point_count();

    if params.enhanced_occupancy_map.is_some() {
        unimplemented!("enhancedOccupanyMap not implemented")
    }
    trace!(
        "frame {}, tile {}, regularPoints+eomPoints {}",
        frame_index,
        tile_index,
        reconstruct.point_count()
    );

    if params.use_additional_points_patch {
        unimplemented!("additionalPointsPatch not implemented")
    }

    if params.geometry_smoothing.is_some() && params.pbf.is_none() {
        unimplemented!("geometrySmoothing && !pbf not implemented")
    }

    for i in 0..attribute_count as usize {
        reconstruct = color_point_cloud(
            reconstruct,
            tile,
            params,
            &atlas.attr_frames[i],
            &point_to_pixel,
            context.get_vps().unwrap().multiple_map_streams_present_flag,
        );
    }

    Some((reconstruct, partition))
}

fn generate_points(
    params: &GeneratePointCloudParams,
    tile: &TileContext,
    geo_video: &VideoGeometry,
    video_frame_index: usize,
    patch_index: usize,
    (u, v): (usize, usize),
    (x, y): (usize, usize),
) -> Vec<Point3D> {
    let patch = &tile.patches[patch_index];
    let frame0 = geo_video.get(video_frame_index).unwrap();
    let mut created_points = vec![];
    let point0 = if params.pbf.is_some() {
        unimplemented!("pbfEnableFlag is not implemented yet");
    } else {
        // DIFF(21Jan23): In this function, we divide the depth by 4 because libavcodec decodes the geometry video as 10-bit (YUV420p10LE), but the reference implementation decodes the geometry video as 8-bit
        // Note that it should still stored as Video<u16> (yeah reference impl is strange i know). We cannot change to Video<u8> otherwise the metadata (indices etc) stored in the atlas will point to out-of-bounds data.
        patch.generate_point(u, v, frame0.get(0, x, y) / 4)
    };
    created_points.push(point0);
    if params.single_map_pixel_interleaving {
        unimplemented!("singleMapPixelInterleaving is not implemented yet");
    } else if params.point_local_reconstruction.is_some() {
        unimplemented!("pointLocalReconstruction is not implemented yet");
    } else if params.map_count_minus1 > 0 {
        let frame1 = if params.multiple_streams {
            unimplemented!("multipleStreams is not implemented yet");
        } else {
            geo_video.get(video_frame_index + 1).unwrap()
        };
        // NOTE: see above for the reason why we divide the depth by 4
        let d1 = frame1.get(0, x, y) / 4_u16;
        let point1 = if params.absolute_d1 {
            patch.generate_point(u, v, d1)
        } else if patch.projection_mode == 0 {
            let mut point1 = point0;
            point1[patch.axes.0 as usize] += d1;
            point1
        } else {
            let mut point1 = point0;
            point1[patch.axes.0 as usize] -= d1;
            point1
        };
        created_points.push(point1);
    }

    created_points
}

/// Color a point set with the attribute.
fn color_point_cloud(
    mut reconstruct: PointSet3,
    tile: &TileContext,
    params: &GeneratePointCloudParams,
    video: &VideoAttribute,
    point_to_pixel: &[Vector3<usize>],
    multiple_streams: bool,
) -> PointSet3 {
    trace!("colorPointCloud start");
    if reconstruct.point_count() == 0 {
        return reconstruct;
    }

    if reconstruct.colors.len() < reconstruct.point_count() {
        reconstruct
            .colors
            .resize_with(reconstruct.point_count(), || Vector3 { x: 0, y: 0, z: 0 });
    }
    assert!(reconstruct.colors.len() == reconstruct.positions.len());

    let _frame0 = video.get(0).unwrap();
    let _frame1 = video.get(1).unwrap();
    let map_count = params.map_count_minus1 as usize + 1;

    let point_count = if tile.use_raw_points_separate_video {
        tile.total_number_of_regular_points
    } else {
        tile.total_number_of_regular_points
            + tile._total_number_of_eom_points
            + tile._total_number_of_raw_points
    };

    trace!("pointcount = {}", point_count);
    trace!(
        "reconstruct.getPointCount() = {}",
        reconstruct.point_count()
    );
    trace!("point_to_pixel size = {}", point_to_pixel.len());
    trace!(
        "pointLocalReconstruction = {:?}",
        params.point_local_reconstruction
    );
    trace!(
        "singleMapPixelInterleaving = {}",
        params.single_map_pixel_interleaving
    );
    trace!(
        "enhancedOccupancyMapCode = {:?}",
        params.enhanced_occupancy_map
    );
    trace!("multipleStreams = {}", multiple_streams);
    let shift = if params.multiple_streams {
        tile.frame_index
    } else {
        tile.frame_index * map_count
    };

    for (i, location) in point_to_pixel.iter().enumerate() {
        let Vector3 { x, y, z } = location.add(Vector3::new(
            tile.left_top_in_frame.0,
            tile.left_top_in_frame.1,
            0,
        ));
        if params.single_map_pixel_interleaving {
            unimplemented!("singleMapPixelInterleaving is not implemented yet");
        } else if multiple_streams {
            unimplemented!("multipleStreams is not implemented yet");
        } else if z < map_count {
            let frame = video.get(z + shift).unwrap();
            for c in 0..3 {
                reconstruct.colors16bit[i][c] = frame.get(c, x, y);
            }
        } else {
            unimplemented!("color transfer weight? is not implemented yet");
        }
    }

    // TODO: Implement color weight transfer here
    // if (target.getPointCount() > 0) { ... }

    if tile.use_raw_points_separate_video {
        unimplemented!("rawPointsSeparateVideo is not implemented yet");
    }

    if tile.use_raw_points_separate_video {
        unimplemented!("rawPointsSeparateVideo is not implemented yet");
    }

    reconstruct
}

/// Convert YUV10bit to RGB8bit
///
/// Modified from https://softpixel.com/~cwright/programming/colorspace/yuv/
fn convert_yuv16_to_rgb8(color16: &Vector3<u16>) -> Vector3<u8> {
    // yuv16 to rgb8
    let clamp = |x: f64| -> u8 {
        if x < 0. {
            0
        } else if x > 255. {
            255
        } else {
            x as u8
        }
    };

    let offset = 512.;
    let scale = 1023.;

    let Vector3 { x: y, y: u, z: v } = *color16;
    let y = y as f64;
    let u = u as f64;
    let v = v as f64;
    let r = y + 1.4075 * (v - offset);
    let g = y - 0.3455 * (u - offset) - (0.7169 * (v - offset));
    let b = y + 1.7790 * (u - offset);
    let r = clamp((r / scale * 255.).floor());
    let g = clamp((g / scale * 255.).floor());
    let b = clamp((b / scale * 255.).floor());
    Vector3 { x: r, y: g, z: b }
}
