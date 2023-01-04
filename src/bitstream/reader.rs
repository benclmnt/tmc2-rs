use crate::{common::context::Context, decoder::PatchOrientation};

// Copied from source/lib/PccLibBitstreamReader/source/PCCBitstreamReader.cpp
use super::{Bitstream, VideoBitstream};
use log::{debug, info, trace};
use num_enum::FromPrimitive;
use std::collections::VecDeque;
use std::rc::Rc;

#[derive(Default)]
struct ReadableAuxData<'a> {
    syntax: Option<&'a mut Context>,
    aux: Option<[usize; 2]>,
}

trait Readable {
    fn from_bitstream(bitstream: &Bitstream, aux_data: ReadableAuxData) -> Self;
}

#[derive(Default)]
struct V3CUnit {
    unit_type: V3CUnitType,
    size: usize,
    bitstream: Bitstream,
}

impl V3CUnit {
    fn peek_type(&self) -> V3CUnitType {
        V3CUnitType::from(self.bitstream.peek(5) as u8)
    }

    /// Originally PCCBitstreamReader::v3cUnit
    fn decode(&self, syntax: &mut Context) -> V3CUnitType {
        // self.bitstream.reset();

        let unit_type = self.decode_header(syntax);
        assert!(self.unit_type == unit_type);
        self.decode_payload(syntax);
        // syntax.bitstream_stat.v3c_unit_type = self.bitstream.size() - position;
        unit_type
    }

    /// Originally PCCBitstreamReader::v3cUnitHeader
    /// V3C Unit Header is 4-byte wide
    fn decode_header(&self, syntax: &mut Context) -> V3CUnitType {
        let v3c_unit_type = V3CUnitType::from(self.bitstream.read(5) as u8); // u(5)

        let mut v3c_unit_header = syntax
            .get_v3c_unit_header(&v3c_unit_type)
            .unwrap_or_default();

        if v3c_unit_type != V3CUnitType::V3cParameterSet {
            v3c_unit_header.sequence_parameter_set_id = self.bitstream.read(4) as u8; // u(4)

            // syntax.active_vps = v3c_unit_header.sequence_parameter_set_id;
            v3c_unit_header.atlas_id = self.bitstream.read(6) as u8; // u(6)
            assert!(
                v3c_unit_header.atlas_id == 0,
                "V3C only have a single atlas"
            );
        }

        match v3c_unit_type {
            V3CUnitType::AttributeVideoData => {
                v3c_unit_header.attribute_index = self.bitstream.read(7) as u8; // u(7)
                v3c_unit_header.attribute_dimension_index = self.bitstream.read(5) as u8; // u(5)
                v3c_unit_header.map_index = self.bitstream.read(4) as u8; // u(4)
                v3c_unit_header.auxiliary_video_flag = self.bitstream.read(1) != 0;
            }
            V3CUnitType::GeometryVideoData => {
                v3c_unit_header.map_index = self.bitstream.read(4) as u8; // u(4)
                v3c_unit_header.auxiliary_video_flag = self.bitstream.read(1) != 0; // u(1)
                self.bitstream.read(12);
            }
            V3CUnitType::OccupancyVideoData | V3CUnitType::AtlasData => {
                self.bitstream.read(17);
            }
            _ => {
                self.bitstream.read(27);
            }
        }

        assert!(
            !v3c_unit_header.auxiliary_video_flag,
            "Auxiliary video not implemented"
        );
        syntax.set_v3c_unit_header(v3c_unit_type, v3c_unit_header);
        v3c_unit_type
    }

    /// Originally PCCBitstreamReader::v3cUnitPayload + videoSubstream
    fn decode_payload(&self, syntax: &mut Context) {
        match self.unit_type {
            V3CUnitType::V3cParameterSet => {
                let vps =
                    V3CParameterSet::from_bitstream(&self.bitstream, ReadableAuxData::default());
                syntax.allocate_atlas_hls(vps.atlas_count_minus1 as usize + 1);
                syntax.add_v3c_parameter_set(vps);
            }
            V3CUnitType::AtlasData => {
                SampleStreamNalUnit::from_bitstream(syntax, &self.bitstream);
            }
            V3CUnitType::OccupancyVideoData => {
                syntax.add_video_bitstream(VideoBitstream::from_bitstream(
                    &self.bitstream,
                    self.size - 4,
                    super::VideoType::Occupancy,
                ));
            }
            V3CUnitType::GeometryVideoData => {
                let vuh = syntax.get_v3c_unit_header(&self.unit_type).unwrap();
                if vuh.auxiliary_video_flag {
                    // syntax.add_video_bitstream(VideoBitstream::from_bitstream(
                    //     &self.bitstream,
                    //     self.size - 4,
                    //     super::VideoType::GeometryRaw,
                    // ));
                    unimplemented!("Auxiliary video for GVD not implemented")
                } else {
                    let vps = syntax.get_vps().unwrap();
                    if vps.map_count_minus1 > 0 && vps.multiple_map_streams_present_flag {
                        unimplemented!("Multiple map streams for GVD not implemented");
                    }
                    syntax.add_video_bitstream(VideoBitstream::from_bitstream(
                        &self.bitstream,
                        self.size - 4,
                        super::VideoType::Geometry,
                    ));
                }
            }
            V3CUnitType::AttributeVideoData => {
                let vuh = syntax.get_v3c_unit_header(&self.unit_type).unwrap();
                let vps = syntax.get_vps().unwrap();
                if vps.attribute_information.attribute_count == 0 {
                    return;
                }
                if vuh.auxiliary_video_flag {
                    unimplemented!("Auxiliary video for AVD not implemented");
                } else if vps.map_count_minus1 > 0 && vps.multiple_map_streams_present_flag {
                    unimplemented!("Multiple map streams for AVD not implemented");
                } else {
                    assert_eq!(
                        vuh.attribute_dimension_index, 0,
                        "attribute_dimension_index > 0 for AVD not implemented"
                    );
                    syntax.add_video_bitstream(VideoBitstream::from_bitstream(
                        &self.bitstream,
                        self.size - 4,
                        super::VideoType::Attribute,
                    ));
                }
            }
        }
    }
}

#[derive(Default, Copy, Clone, Debug)]
pub(crate) struct V3CUnitHeader {
    /// value of vps_v3c_parameter_set_id for the active v3c VPS
    sequence_parameter_set_id: u8,
    /// ID of the atlas that corresponds to the current V3C unit
    atlas_id: u8,
    /// index of the attribute data carried in AVD         
    attribute_index: u8,
    /// index of the attribute dimension group carried in AVD
    attribute_dimension_index: u8,
    /// when present indicates the map index of the current geometry or attribute stream.
    /// when absent, the map index is derived from the type of the sub-bitstream
    map_index: u8,
    /// indicates if the associated GVD / AVD is a RAW and/or EOM coded points video-only sub-bitstream
    auxiliary_video_flag: bool,
}

#[derive(Default, PartialEq, Eq, Hash, Copy, Clone, Debug, FromPrimitive)]
#[repr(u8)]
pub(crate) enum V3CUnitType {
    #[default]
    /// Sequence parameter set
    V3cParameterSet,
    /// Patch Data Group / Atlas information
    AtlasData,
    OccupancyVideoData,
    GeometryVideoData,
    AttributeVideoData,
    // Packed Video Data
    // PackedVideoData,
}

#[derive(Default)]
pub struct V3CParameterSet {
    pub(crate) profile_tier_level: ProfileTierLevel,
    pub(crate) v3c_parameter_set_id: u8,
    /// (12Dec22) This is always 0 up to reference impl v18.
    /// There is a comment in the original code (PCCEncoderParameters.cpp) saying V3C only have a single atlas
    ///
    /// DIFF: Making use of the fact above, we set all the following fields to its non-Vec equivalent
    pub(crate) atlas_count_minus1: u8,
    atlas_id: u8, // originally Vec<u16>
    pub(crate) frame_width: u16,
    frame_height: u16,
    pub(crate) map_count_minus1: u8,
    pub(crate) multiple_map_streams_present_flag: bool,
    pub(crate) map_absolute_coding_enable_flag: Vec<bool>,
    map_predictor_index_diff: Vec<bool>,
    pub(crate) auxiliary_video_present_flag: bool,
    pub(crate) occupancy_video_present_flag: bool,
    pub(crate) geometry_video_present_flag: bool,
    pub(crate) attribute_video_present_flag: bool,
    pub(crate) occupancy_information: OccupancyInformation,
    pub(crate) geometry_information: GeometryInformation,
    pub(crate) attribute_information: AttributeInformation,
    extension_present_flag: bool,
    // extension_8bits: u8,
    // extension_length_minus1: usize,
    // extension_data_byte: Vec<u8>,
    // vps_vpcc_extension: VpcVpccExtension,
}

impl V3CParameterSet {
    /// Allocate atlas data
    /// All implementation is commented because atlas_count_minus1 is always 0.
    #[inline]
    fn allocate_atlas(&mut self) {
        // self.atlas_id
        //     .resize(self.atlas_count_minus1 as usize + 1, 0);
        // self.frame_width
        //     .resize(self.atlas_count_minus1 as usize + 1, 1);
        // self.frame_height
        //     .resize(self.atlas_count_minus1 as usize + 1, 0);
        // self.map_count_minus1
        //     .resize(self.atlas_count_minus1 as usize + 1, 0);
        // self.multiple_map_streams_present_flag
        //     .resize(self.atlas_count_minus1 as usize + 1, false);
        // self.map_absolute_coding_enable_flag
        //     .resize_with(self.atlas_count_minus1 as usize + 1, Vec::new);
        // self.map_predictor_index_diff
        //     .resize_with(self.atlas_count_minus1 as usize + 1, Vec::new);
        // self.auxiliary_video_present_flag
        //     .resize(self.atlas_count_minus1 as usize + 1, false);
        // self.occupancy_video_present_flag
        //     .resize(self.atlas_count_minus1 as usize + 1, false);
        // self.geometry_video_present_flag
        //     .resize(self.atlas_count_minus1 as usize + 1, false);
        // self.attribute_video_present_flag
        //     .resize(self.atlas_count_minus1 as usize + 1, false);
        // self.occupancy_information.resize_with(
        //     self.atlas_count_minus1 as usize + 1,
        //     OccupancyInformation::default,
        // );
        // self.geometry_information.resize_with(
        //     self.atlas_count_minus1 as usize + 1,
        //     GeometryInformation::default,
        // );
        // self.attribute_information.resize_with(
        //     self.atlas_count_minus1 as usize + 1,
        //     AttributeInformation::default,
        // );
    }

    /// allocate maps at index `index`
    fn allocate_map(&mut self) {
        self.map_absolute_coding_enable_flag
            .resize(self.map_count_minus1 as usize + 1, true);
        self.map_predictor_index_diff
            .resize(self.map_count_minus1 as usize + 1, false);
    }
}

impl Readable for V3CParameterSet {
    fn from_bitstream(bitstream: &Bitstream, _: ReadableAuxData) -> Self {
        let mut sps = V3CParameterSet {
            profile_tier_level: ProfileTierLevel::from_bitstream(bitstream),
            ..Default::default()
        };
        debug!("[ptl] {:?}", &sps.profile_tier_level);
        sps.v3c_parameter_set_id = bitstream.read(4) as u8; // u(4)
        bitstream.read(8); // u(8)
        sps.atlas_count_minus1 = bitstream.read(6) as u8; // u(6)
        assert!(
            sps.atlas_count_minus1 == 0,
            "V3C only have a single atlas. This is the case up to mpeg-pcc-tmc2 v18."
        );
        sps.allocate_atlas();

        sps.atlas_id = bitstream.read(6) as u8; // u(6)
        sps.frame_width = bitstream.read_uvlc() as u16; // ue(v)
        sps.frame_height = bitstream.read_uvlc() as u16; // ue(v)
        sps.map_count_minus1 = bitstream.read(4) as u8; // u(4)
        sps.allocate_map();
        if sps.map_count_minus1 > 0 {
            sps.multiple_map_streams_present_flag = bitstream.read(1) != 0;
            // u(1)
            assert!(
                !sps.multiple_map_streams_present_flag,
                "V3C only have a single map. This is the case up to mpeg-pcc-tmc2 v18."
            )
        }
        sps.map_absolute_coding_enable_flag.push(true);
        for k in 1..=sps.map_count_minus1 as usize {
            if sps.multiple_map_streams_present_flag {
                sps.map_absolute_coding_enable_flag[k] = bitstream.read(1) != 0;
                // u(1)
            }

            if !sps.map_absolute_coding_enable_flag[k] {
                sps.map_predictor_index_diff[k] = bitstream.read_uvlc() != 0;
                // ue(v)
            }
        }
        sps.auxiliary_video_present_flag = bitstream.read(1) != 0; // u(1)
        sps.occupancy_video_present_flag = bitstream.read(1) != 0; // u(1)
        sps.geometry_video_present_flag = bitstream.read(1) != 0; // u(1)
        sps.attribute_video_present_flag = bitstream.read(1) != 0;
        // u(1)
        trace!(
            "Atlas *_video_present flags: aux: {:?} occupancy: {:?} geom: {:?} attr: {:?}",
            &sps.auxiliary_video_present_flag,
            &sps.occupancy_video_present_flag,
            &sps.geometry_video_present_flag,
            &sps.attribute_video_present_flag
        );

        if sps.occupancy_video_present_flag {
            sps.occupancy_information = OccupancyInformation::from_bitstream(bitstream);
        }
        if sps.geometry_video_present_flag {
            sps.geometry_information =
                GeometryInformation::from_bitstream(bitstream, sps.auxiliary_video_present_flag);
        }
        if sps.attribute_video_present_flag {
            sps.attribute_information = AttributeInformation::from_bitstream(
                bitstream,
                sps.auxiliary_video_present_flag,
                sps.map_count_minus1,
            );
        }

        sps.extension_present_flag = bitstream.read(1) != 0; // u(1)
        if sps.extension_present_flag {
            unimplemented!("extension_present_flag");
            // sps.extension_8bits = bitstream.read(8) as u8; // u(8)
        }
        // if sps.extension_8bits > 0 {
        //     sps.extension_length_minus1 = bitstream.read_uvlc() as usize;
        //     sps.extension_data_byte = vec![0; sps.extension_length_minus1 + 1];
        //     for i in 0..sps.extension_length_minus1 + 1 {
        //         sps.extension_data_byte[i] = bitstream.read(8) as u8;
        //         // u(8)
        //     }
        // }
        bitstream.byte_align();
        sps
    }
}

/// 8.3.4.3 Occupancy information Set Syntax
#[derive(Debug)]
pub(crate) struct OccupancyInformation {
    pub(crate) occupancy_codec_id: u8,
    pub(crate) occupancy_lossy_compression_threshold: u8,
    pub(crate) occupancy_2d_bitdepth_minus1: u8,
    pub(crate) occupancy_msb_align_flag: bool,
}

impl Default for OccupancyInformation {
    fn default() -> Self {
        Self {
            occupancy_codec_id: 0,
            occupancy_lossy_compression_threshold: 0,
            occupancy_2d_bitdepth_minus1: 10,
            occupancy_msb_align_flag: false,
        }
    }
}

impl OccupancyInformation {
    fn from_bitstream(bitstream: &Bitstream) -> Self {
        Self {
            occupancy_codec_id: bitstream.read(8) as u8,
            occupancy_lossy_compression_threshold: bitstream.read(8) as u8,
            occupancy_2d_bitdepth_minus1: bitstream.read(5) as u8,
            occupancy_msb_align_flag: bitstream.read(1) != 0,
        }
    }
}

/// 8.3.4.4 Geometry information Syntax
#[derive(Debug)]
pub(crate) struct GeometryInformation {
    pub(crate) geometry_codec_id: u8,
    auxiliary_geometry_codec_id: u8,
    pub(crate) geometry_2d_bitdepth_minus1: u8,
    pub(crate) geometry_3d_coordinates_bitdepth_minus1: u8,
    pub(crate) geometry_msb_align_flag: bool,
}

impl Default for GeometryInformation {
    fn default() -> Self {
        Self {
            geometry_codec_id: 0,
            auxiliary_geometry_codec_id: 0,
            geometry_2d_bitdepth_minus1: 10,
            geometry_3d_coordinates_bitdepth_minus1: 9,
            geometry_msb_align_flag: false,
        }
    }
}

impl GeometryInformation {
    fn from_bitstream(bitstream: &Bitstream, is_auxiliary_video_present: bool) -> Self {
        Self {
            geometry_codec_id: bitstream.read(8) as u8, // u(8)
            geometry_2d_bitdepth_minus1: bitstream.read(5) as u8, // u(5)
            geometry_msb_align_flag: bitstream.read(1) != 0, // u(1)
            geometry_3d_coordinates_bitdepth_minus1: bitstream.read(5) as u8, // u(5)
            auxiliary_geometry_codec_id: if is_auxiliary_video_present {
                bitstream.read(8) as u8 // u(8)
            } else {
                0
            },
        }
    }
}

/// 8.3.4.5 Attribute Information Syntax
/// TODO: can be optimized since we only have 1 attribute
#[derive(Debug, Default)]
pub(crate) struct AttributeInformation {
    pub(crate) attribute_count: u8,
    pub(crate) attribute_type_id: Vec<u8>,
    pub(crate) attribute_codec_id: Vec<u8>,
    pub(crate) auxiliary_attribute_codec_id: Vec<u8>,
    pub(crate) attribute_map_absolute_coding_persistence_flag: Vec<bool>,
    pub(crate) attribute_dimension_minus1: Vec<u8>,
    pub(crate) attribute_dimension_partitions_minus1: Vec<u8>,
    pub(crate) attribute_partition_channels_minus1: Vec<Vec<u8>>,
    pub(crate) attribute_2d_bitdepth_minus1: Vec<u8>,
    pub(crate) attribute_msb_align_flag: Vec<bool>,
}

impl AttributeInformation {
    fn new(attribute_count: usize) -> Self {
        Self {
            attribute_count: attribute_count as u8,
            attribute_type_id: vec![0; attribute_count],
            attribute_codec_id: vec![0; attribute_count],
            auxiliary_attribute_codec_id: vec![0; attribute_count],
            attribute_map_absolute_coding_persistence_flag: vec![false; attribute_count],
            attribute_dimension_minus1: vec![0; attribute_count],
            attribute_dimension_partitions_minus1: vec![0; attribute_count],
            attribute_partition_channels_minus1: vec![Vec::new(); attribute_count],
            attribute_2d_bitdepth_minus1: vec![0; attribute_count],
            attribute_msb_align_flag: vec![false; attribute_count],
        }
    }

    fn from_bitstream(
        bitstream: &Bitstream,
        is_auxiliary_video_present: bool,
        map_count_minus1: u8,
    ) -> Self {
        let attribute_count = bitstream.read(7) as usize;
        let mut ai = AttributeInformation::new(attribute_count);
        for i in 0..attribute_count {
            ai.attribute_type_id[i] = bitstream.read(4) as u8; // u(4)
            ai.attribute_codec_id[i] = bitstream.read(8) as u8; // u(8)
            if is_auxiliary_video_present {
                ai.auxiliary_attribute_codec_id[i] = bitstream.read(8) as u8 // u(8)
            }
            ai.attribute_map_absolute_coding_persistence_flag[i] = true;
            if map_count_minus1 > 0 {
                ai.attribute_map_absolute_coding_persistence_flag[i] = bitstream.read(1) != 0;
            }
            ai.attribute_dimension_minus1[i] = bitstream.read(6) as u8; // u(8)
            if ai.attribute_dimension_minus1[i] > 0 {
                ai.attribute_dimension_partitions_minus1[i] = bitstream.read(6) as u8; // u(6)
                let mut remaining_dimensions = ai.attribute_dimension_minus1[i];
                let k = ai.attribute_dimension_partitions_minus1[i];
                for j in 0..k as usize {
                    let partition_channels = if k - (j as u8) == remaining_dimensions {
                        0
                    } else {
                        bitstream.read_uvlc() as u8
                    };
                    ai.attribute_partition_channels_minus1[i].push(partition_channels);
                    remaining_dimensions -= partition_channels
                }
                ai.attribute_partition_channels_minus1[i].push(remaining_dimensions);
            }
            ai.attribute_2d_bitdepth_minus1[i] = bitstream.read(5) as u8;
            ai.attribute_msb_align_flag[i] = bitstream.read(1) != 0;
        }
        ai
    }
}

/// 8.3.4.2 Profile, Tier, and Level syntax
///
/// Profiles, tiers, and levels specify restrictions on the bitstreams and hence limits on the capabilities needed to decode the bitstreams.
/// Profiles, tiers, and levels may also be used to indicate interoperability points between individual decoder implementations.
///
/// Each profile specifies a subset of algorithmic features and limits that shall be supported by all decoders conforming to that profile
/// Each level of a tier specifies a set of limits on the values that may be taken by the syntax elements of the bitstream. Generally corresponds to a decoder's processing load and memory capability
///
/// A V-PCC profile consists of a combination of up to three profile components, identified from the syntax elements
/// ptl_profile_codec_group_idc, ptl_profile_pcc_toolset_idc, and ptl_profile_reconstruction_idc
///     - The first conformance point, point A, covers the decoded attributes, geometry, and occupancy video sub-bitstreams, plus the decoded atlas sub-bitstream.
///       It also covers the derived block to patch map information, but doesn't cover the point cloud reconstruction process.
///       ptl_profile_codec_group_idc and ptl_profile_pcc_toolset_idc information describe conformance point A.
///     - The second conformance point, point B, covers the point cloud reconstruction process.
///       ptl_profile_reconstruction_idc information describes conformance point B.
///
#[derive(Debug, Default)]
pub(crate) struct ProfileTierLevel {
    /// Main (0)
    tier_flag: bool,

    /// AVC Progressive High, HEVC Main 10, HEVC444
    pub(crate) profile_codec_group_idc: u8,
    /// Indicates the use of V-PCC specific tools
    /// Basic (0) or Extended (1)
    profile_toolset_idc: u8,
    /// Rec0 (0) / Rec1 (1) / Rec Unconstrained (2)
    /// Rec0: ignore pixel deinterleaving, PLR, EOM, duplicate point removal, RAW patches, smoothing
    /// Rec1: reconstruct pixel deinterleaving, PLR, EOM, duplicate point removal, RAW patches, smoothing
    profile_reconstruction_idc: u8,

    /// Level 1.0 / 2.0
    level_idc: u8,
    // extended_sub_profile_flag: bool,
    // num_sub_profiles: u8,
    // sub_profile_idc: Vec<u8>,
    // tool_constraints_present_flag: bool,
    // profile_toolset_constraints_information: ProfileToolsetConstraintsInformation,
}

impl ProfileTierLevel {
    fn from_bitstream(bitstream: &Bitstream) -> Self {
        let mut ptl = ProfileTierLevel {
            tier_flag: bitstream.read(1) != 0,                   // u(1)
            profile_codec_group_idc: bitstream.read(7) as u8,    // u(7)
            profile_toolset_idc: bitstream.read(8) as u8,        // u(8)
            profile_reconstruction_idc: bitstream.read(8) as u8, // u(8)
            ..Default::default()
        };
        // reserved 32 zero bits
        bitstream.move_to_next_byte();
        bitstream.move_to_next_byte();
        bitstream.move_to_next_byte();
        bitstream.move_to_next_byte();

        ptl.level_idc = bitstream.read(8) as u8; // u(8)

        // ptl.num_sub_profiles =
        assert!(
            bitstream.read(6) as u8 == 0,
            "(9Dec22) ptl subprofiles not supported"
        ); // u(6)
           // ptl.extended_sub_profile_flag =
        let _ = bitstream.read(1) != 0; // u(1)

        // ptl.sub_profile_idc
        //     .reserve_exact(ptl.num_sub_profiles as usize);
        // for _ in 0..ptl.num_sub_profiles {
        //     assert!(!ptl.extended_sub_profile_flag);
        //     ptl.sub_profile_idc.push(bitstream.read(32) as u8); // u(8)
        // }

        // ptl.tool_constraints_present_flag =
        assert!(
            bitstream.read(1) == 0,
            "(9Dec22) ptl toolset constraints information not supported"
        ); // u(1)

        // if ptl.tool_constraints_present_flag {
        //     unimplemented!("profile toolset constraints information unimplemented")
        //     // ptl.profile_toolset_constraints_information =
        //     //     ProfileToolsetConstraintsInformation::from_bitstream(bitstream);
        // }
        ptl
    }
}

#[derive(Debug, Default)]
struct ProfileToolsetConstraintsInformation {
    one_frame_only_flag: bool,
    eom_constraint_flag: bool,
    max_map_count_minus1: u8,
    max_atlas_count_minus1: u8,
    multiple_map_streams_constraint_flag: bool,
    plr_constraint_flag: bool,
    attribute_max_dimension_minus1: u8,
    attribute_max_dimension_partitions_minus1: u8,
    no_eight_orientations_constraint_flag: bool,
    no_45deg_projection_patch_constraint_flag: bool,
    num_reserved_constraint_bytes: u8,
    reserved_constraint_bytes: Vec<u8>,
}

impl ProfileToolsetConstraintsInformation {
    pub fn from_bitstream(bitstream: &Bitstream) -> Self {
        let mut ptci = ProfileToolsetConstraintsInformation {
            one_frame_only_flag: bitstream.read(1) != 0,     // u(1)
            eom_constraint_flag: bitstream.read(1) != 0,     // u(1)
            max_map_count_minus1: bitstream.read(4) as u8,   // u(4)
            max_atlas_count_minus1: bitstream.read(4) as u8, // u(4)
            multiple_map_streams_constraint_flag: bitstream.read(1) != 0, // u(1)
            plr_constraint_flag: bitstream.read(1) != 0,     // u(1)
            attribute_max_dimension_minus1: bitstream.read(6) as u8, // u(6)
            attribute_max_dimension_partitions_minus1: bitstream.read(6) as u8, // u(6)

            no_eight_orientations_constraint_flag: bitstream.read(1) != 0, // u(1)
            no_45deg_projection_patch_constraint_flag: bitstream.read(1) != 0, // u(1)
            ..Default::default()
        };
        bitstream.byte_align();

        ptci.num_reserved_constraint_bytes = bitstream.read(8) as u8; // u(8)
        ptci.reserved_constraint_bytes
            .reserve_exact(ptci.num_reserved_constraint_bytes as usize);
        for _ in 0..ptci.num_reserved_constraint_bytes {
            ptci.reserved_constraint_bytes.push(bitstream.read(8) as u8); // u(8)
        }
        ptci
    }
}

pub struct SampleStreamV3CUnit {
    units: VecDeque<V3CUnit>,
    pub(super) ssvh_unit_size_precision_bytes_minus1: u8,
}

impl SampleStreamV3CUnit {
    /// Read v3c units from bitstream into ssvu
    /// Originally: PCCBitstreamReader::read
    pub fn from_bitstream(bitstream: &Bitstream) -> (Self, usize) {
        let mut ssvu = Self {
            units: VecDeque::new(),
            ssvh_unit_size_precision_bytes_minus1: SampleStreamV3CUnit::read_header(bitstream),
        };
        let mut header_size = 1;

        while bitstream.more_data() {
            let v3c_unit = SampleStreamV3CUnit::read_v3c_unit(
                bitstream,
                ssvu.ssvh_unit_size_precision_bytes_minus1 + 1,
            );
            header_size += ssvu.ssvh_unit_size_precision_bytes_minus1 as usize + 1;
            ssvu.units.push_back(v3c_unit);
        }
        (ssvu, header_size)
    }

    /// Gets the precision bytes minus 1 for ssvh unit size from the bitstream
    #[inline]
    fn read_header(bitstream: &Bitstream) -> u8 {
        let precision_bytes_minus1 = bitstream.read(3) as u8;
        bitstream.read(5); // padding i think
        precision_bytes_minus1
    }

    /// Reads a V3CUnit from the bitstream
    fn read_v3c_unit(bitstream: &Bitstream, precision_bytes: u8) -> V3CUnit {
        let mut v3c_unit = V3CUnit {
            size: bitstream.read(8 * precision_bytes) as usize,
            ..Default::default()
        };
        v3c_unit
            .bitstream
            .copy_from(bitstream, bitstream.bytes(), v3c_unit.size);
        // DIFF: reset v3c unit here for decode function below
        v3c_unit.bitstream.reset();
        let v3c_unit_type8 = v3c_unit.bitstream.data[0] >> 3;
        let v3c_unit_type: V3CUnitType = v3c_unit_type8.into();
        v3c_unit.unit_type = v3c_unit_type;
        info!(
            "[v3c_unit] size: {}, type: {:?}",
            v3c_unit.size, v3c_unit.unit_type
        );
        v3c_unit
    }

    /// Originally PCCBitstreamReader::decode
    pub fn decode(&mut self, syntax: &mut Context) {
        let mut num_vps = 0; // counter for atlas information

        // let bitstream_stat = &mut syntax.bitstream_stat;
        // bitstream_stat.new_gof();

        while self.get_v3c_unit_count() > 0 {
            let unit = &mut self.front();
            let v3c_unit_type = unit.peek_type();
            if v3c_unit_type == V3CUnitType::V3cParameterSet {
                num_vps += 1;
                if num_vps > 1 {
                    // remove the bits counted for the last VPS
                    // let v3c_unit_size = unit.bitstream.size();
                    // TODO[stat]: let stat_size = bitstream_stat.overwrite_v3c_unit_size

                    // Each V3C Unit only contains 1 VPS. In the reference implementation, we trackback.
                    // However, this causes `unit.decode(syntax)` line above to be evaluated twice.
                    // DIFF: Instead, we could peek first to determine if it is a VPS

                    // Note: the usage of endOfGop is changed to break
                    break;
                }
            }
            unit.decode(syntax);
            self.pop_front();
        }
    }

    fn front(&self) -> &V3CUnit {
        self.units.front().unwrap()
    }

    fn pop_front(&mut self) {
        self.units.pop_front();
    }

    pub fn get_v3c_unit_count(&self) -> usize {
        self.units.len()
    }
}

struct SampleStreamNalUnit {
    units: VecDeque<NalUnit>,
    size_precision_bytes_minus1: u8,
}

impl SampleStreamNalUnit {
    pub fn from_bitstream(syntax: &mut Context, bitstream: &Bitstream) -> Self {
        let mut ssnu = Self {
            units: VecDeque::new(),
            size_precision_bytes_minus1: SampleStreamNalUnit::read_header(bitstream),
        };
        let mut prefix_sei: Rc<Option<SeiRbsp>> = Rc::new(None);
        while bitstream.more_data() {
            // Originally: PCCBitstreamReader::sampleStreamNalUnit
            let nalu_size = bitstream.read(8 * (ssnu.size_precision_bytes_minus1 + 1)) as usize;
            let nalu = NalUnit::from_bitstream(syntax, bitstream, nalu_size, &mut prefix_sei);
            debug!(
                "[nalu] size: {}, precision: {}, type: {:?}, bitsWritten: {}",
                nalu.size,
                ssnu.size_precision_bytes_minus1 + 1,
                nalu.unit_type,
                bitstream.position_in_bytes()
            );
            ssnu.units.push_back(nalu);
        }
        ssnu
    }

    /// Originally: PCCBitstreamReader::sampleStreamNalHeader
    #[inline]
    fn read_header(bitstream: &Bitstream) -> u8 {
        let precision_bytes_minus1 = bitstream.read(3) as u8;
        bitstream.read(5); // padding i think
        precision_bytes_minus1
    }
}

/// 8.3.5 NAL (Network Abstraction Layer) unit syntax
struct NalUnit {
    unit_type: NalUnitType,
    layer_id: u8,
    temporal_id_plus_1: u8,
    size: usize,
    // data: Vec<u8>,
}

impl NalUnit {
    pub fn from_bitstream(
        syntax: &mut Context,
        bitstream: &Bitstream,
        nalu_size: usize,
        prefix_sei: &mut Rc<Option<SeiRbsp>>,
    ) -> Self {
        // The part of PCCBitstreamReader::sampleStreamNalUnit before the call to nalUnitHeader is moved out to the caller.

        // originally PCCBitstreamReader::nalUnitHeader
        bitstream.read(1);
        let nalu = Self {
            size: nalu_size,
            unit_type: NalUnitType::from(bitstream.read(6) as u8), // u(6)
            layer_id: bitstream.read(6) as u8,                     // u(6)
            temporal_id_plus_1: bitstream.read(3) as u8,           // u(3)
                                                                   // data: Vec::with_capacity(nalu_size - 2),
        };

        // continue to PCCBitstreamReader::sampleStreamNalUnit
        match nalu.unit_type {
            NalUnitType::Asps => {
                syntax.add_atlas_sequence_parameter_set(
                    AtlasSequenceParameterSetRbsp::from_bitstream(bitstream),
                );
            }
            NalUnitType::Afps => {
                syntax.add_atlas_frame_parameter_set(AtlasFrameParameterSetRbsp::from_bitstream(
                    syntax, bitstream,
                ));
            }
            NalUnitType::TrailN
            | NalUnitType::TrailR
            | NalUnitType::TsaN
            | NalUnitType::TsaR
            | NalUnitType::StsaN
            | NalUnitType::StsaR
            | NalUnitType::RadlN
            | NalUnitType::RadlR
            | NalUnitType::RaslN
            | NalUnitType::RaslR
            | NalUnitType::SkipN
            | NalUnitType::SkipR
            | NalUnitType::IdrNLp => {
                let mut atl = AtlasTileLayerRbsp::from_bitstream(syntax, bitstream, nalu.unit_type);
                atl.sei = prefix_sei.clone();
                syntax.add_atlas_tile_layer(atl);
            }
            NalUnitType::PrefixESEI | NalUnitType::PrefixNSEI => {
                let sei = SeiRbsp::from_bitstream(bitstream, nalu.unit_type);
                *Rc::get_mut(prefix_sei).unwrap() = Some(sei);
            }
            NalUnitType::SuffixESEI | NalUnitType::SuffixNSEI => {
                unimplemented!("suffixSEI not implemented")
            }
            _ => unreachable!("Unknown NAL unit type"),
        };
        nalu
    }
}

/// ACL (Atlas Coding Layer):
#[derive(Debug, PartialEq, FromPrimitive, Clone, Copy)]
#[repr(u8)]
pub(crate) enum NalUnitType {
    ///  0-1: Coded tile of a non-TSA, non-STSA trailing atlas frame, ACL
    #[default]
    TrailN,
    TrailR,

    ///  2-3: Coded tile of a TSA (Temporal Sub-layer Access) atlas frame, ACL
    TsaN,
    TsaR,

    /// 4-5: Coded tile of a STSA (Step-wise Temporal Sub-layer Access) atlas frame, ACL
    StsaN,
    StsaR,

    /// 6-7: Coded tile of a RADL (Random Access Decodable Leading) atlas frame ACL
    /// RADL coded atlases are NOT used as reference atlases for the decoding process of trailing coded atlases of the same associated IRAP coded atlas.
    /// When present, all RADL coded atlases precede, in decoding order, all trailing coded atlases of the same associated IRAP coded atlas.
    RadlN,
    RadlR,

    /// 8-9: Coded tile of a RASL (Random Access Skipped Leading) atlas frame ACL
    /// All RASL coded atlases are leading coded atlas of an associated BLA or CRA coded atlas.
    /// When the associated IRAP coded atlas has NoOutputBeforeRecoveryFlag equal to 1, the RASL coded atlas is not output and may not be correctly decodable,
    /// as the RASL coded atlas may contain references to coded atlases that are not present in the bitstream.
    /// RASL coded atlases are not used as reference atlases for the decoding process of non-RASL coded atlases.
    /// When present, all RASL coded atlases precede, in decoding order, all trailing coded atlases of the same associated IRAP coded atlas.
    RaslN,
    RaslR,
    /// 10-11: Coded tile of a skipped atlas frame ACL
    SkipN,
    SkipR,

    /// 12,14: Reserved non-IRAP sub-layer non-reference ACL NAL unit types ACL
    /// 13,15: Reserved non-IRAP sub-layer reference ACL NAL unit types ACL
    // RsvAclN12,
    // RsvAclR13,
    // RsvAclN14,
    // RsvAclR15,

    /// 16-18: Coded tile of a BLA (Broken Link Access) atlas frame, ACL
    ///
    /// A BLA coded atlas does NOT use inter prediction in its decoding process, i.e. I-frame
    /// Each BLA coded atlas begins a new CAS (Coded Atlas Sequence) and
    /// it has the same effect on the decoding process as an instantaneous decoding refresh (IDR) coded atlas.
    BlaWLp = 16,
    BlaWRadl,
    BlaNLp,
    /// 19-21: Coded tile of a GBLA (Global Broken Link Access) atlas frame, ACL
    GblaWLp,
    GblaWRadl,
    GblaNLp,

    /// 22-23: Coded tile of a IDR (Instantaneous Decoding Refresh) atlas frame ACL
    /// An IDR Coded atlas does not refer to any atlases other than itself for inter prediction in its decoding process, i.e. I-frame
    IdrWRadl,
    IdrNLp,
    /// 24-25: Coded tile of a GIDR(Global Instantaneous Decoding Refresh) atlas frame, ACL
    GidrWRadl,
    GidrNLp,

    /// 26: Coded tile of a CRA (Clean Random Access) atlas frame, ACL.
    ///
    /// A CRA coded atlas does not use inter prediction in its decoding process, i.e. I-frame,
    /// It could have associated RADL / RASL-coded atlas frames
    Cra,
    /// 27: Coded tile of a GCRA (Global Clean Random Access) atlas frame, ACL
    /// access unit in which the coded atlas (3.34) with nal_layer_id equal to 0 is a GCRA coded atlas (3.63)
    Gcra,

    /// 28-29: Reserved IRAP (Intra Random Access Point) ACL NAL unit types. ACL
    // RsvIrapAcl28,
    // RsvIrapAcl29,
    /// 30-35: Reserved non-IRAP ACL NAL unit types, ACL
    // RsvAcl30,
    // RsvAcl31,
    // RsvAcl32,
    // RsvAcl33,
    // RsvAcl34,
    // RsvAcl35,

    /// 36: Atlas sequence parameter set, non-ACL
    Asps = 36,

    /// 37: Atlas frame parameter set, non-ACL
    Afps,

    /// 38: Access unit delimiter, non-ACL
    Aud,
    /// 39: V3C access unit delimiter, non-ACL
    V3cAud,

    /// 40: End of sequence, non-ACL
    Eos,
    /// 41: End of bitstream, non-ACL        
    Eob,
    /// 42: Filler data, non-ACL         
    Fd,

    /// 43-44: Non-essential supplemental enhancement information, non-ACL         
    PrefixNSEI,
    SuffixNSEI,
    /// 45-46: Essential supplemental enhancement information, non-ACL         
    PrefixESEI,
    SuffixESEI,

    /// 47: Atlas adaptation parameter set, non-ACL
    Aaps,
}

impl NalUnitType {
    fn is_prefix_sei(&self) -> bool {
        matches!(self, NalUnitType::PrefixNSEI | NalUnitType::PrefixESEI)
    }

    fn is_suffix_sei(&self) -> bool {
        matches!(self, NalUnitType::SuffixNSEI | NalUnitType::SuffixESEI)
    }
}

/// Rbsp: Raw bit payload
pub(crate) struct AtlasSequenceParameterSetRbsp {
    atlas_sequence_parameter_set_id: u8,
    pub(crate) frame_width: u16,
    pub(crate) frame_height: u16,
    pub(crate) geometry_2d_bitdepth_minus1: u8,
    pub(crate) geometry_3d_bitdepth_minus1: u8,
    pub(crate) log2_max_atlas_frame_order_cnt_lsb_minus_4: u8,
    max_dec_atlas_frame_buffering_minus1: u8,
    long_term_ref_atlas_frames_flag: bool,
    num_ref_atlas_frame_lists_in_asps: u8,
    pub(crate) ref_list_struct: Vec<RefListStruct>,
    use_eight_orientations_flag: bool,
    pub(crate) extended_projection_enabled_flag: bool,
    max_number_projections_minus1: usize,
    normal_axis_limits_quantization_enabled_flag: bool,
    normal_axis_max_delta_value_enabled_flag: bool,
    patch_precedence_order_flag: bool,
    pub(crate) log2_patch_packing_block_size: u8,
    pub(crate) patch_size_quantizer_present_flag: bool,
    map_count_minus1: u8,
    pub(crate) pixel_deinterleaving_flag: bool,
    // pixel_deinterleaving_map_flag: Vec<bool>,
    pub(crate) eom_patch_enabled_flag: bool,
    pub(crate) eom_fix_bit_count_minus1: u8,
    pub(crate) raw_patch_enabled_flag: bool,
    pub(crate) auxiliary_video_enabled_flag: bool,
    pub(crate) plr_enabled_flag: bool,
    // plr_information: Vec<PLRInformation>,
    vui_parameters_present_flag: bool,
    extension_flag: bool,
    vpcc_extension_flag: bool,
    extension_7bits: u8,
    // vui_parameters: VUIParameters,
    pub(crate) vpcc_extension: AspsVpccExtension,
}

impl Default for AtlasSequenceParameterSetRbsp {
    fn default() -> Self {
        AtlasSequenceParameterSetRbsp {
            atlas_sequence_parameter_set_id: 0,
            frame_width: 0,
            frame_height: 0,
            geometry_2d_bitdepth_minus1: 0,
            geometry_3d_bitdepth_minus1: 0,
            log2_max_atlas_frame_order_cnt_lsb_minus_4: 4,
            max_dec_atlas_frame_buffering_minus1: 0,
            long_term_ref_atlas_frames_flag: false,
            num_ref_atlas_frame_lists_in_asps: 0,
            ref_list_struct: Vec::new(),
            use_eight_orientations_flag: false,
            extended_projection_enabled_flag: false,
            max_number_projections_minus1: 5,
            normal_axis_limits_quantization_enabled_flag: true,
            normal_axis_max_delta_value_enabled_flag: false,
            patch_precedence_order_flag: false,
            log2_patch_packing_block_size: 0,
            patch_size_quantizer_present_flag: false,
            map_count_minus1: 0,
            pixel_deinterleaving_flag: false,
            // pixel_deinterleaving_map_flag: Vec::new(),
            eom_patch_enabled_flag: false,
            eom_fix_bit_count_minus1: 0,
            raw_patch_enabled_flag: false,
            auxiliary_video_enabled_flag: false,
            plr_enabled_flag: false,
            // plr_information: Vec::new(),
            vui_parameters_present_flag: false,
            extension_flag: false,
            vpcc_extension_flag: false,
            extension_7bits: 0,
            // vui_parameters: VUIParameters::default(),
            vpcc_extension: AspsVpccExtension::default(),
        }
    }
}

impl AtlasSequenceParameterSetRbsp {
    fn from_bitstream(bitstream: &Bitstream) -> Self {
        let mut asps = AtlasSequenceParameterSetRbsp::default();
        asps.atlas_sequence_parameter_set_id = bitstream.read_uvlc() as u8;
        asps.frame_width = bitstream.read_uvlc() as u16;
        asps.frame_height = bitstream.read_uvlc() as u16;
        asps.geometry_3d_bitdepth_minus1 = bitstream.read(5) as u8;
        asps.geometry_2d_bitdepth_minus1 = bitstream.read(5) as u8;
        asps.log2_max_atlas_frame_order_cnt_lsb_minus_4 = bitstream.read_uvlc() as u8;
        asps.max_dec_atlas_frame_buffering_minus1 = bitstream.read_uvlc() as u8;
        asps.long_term_ref_atlas_frames_flag = bitstream.read(1) != 0;
        asps.num_ref_atlas_frame_lists_in_asps = bitstream.read_uvlc() as u8;
        asps.ref_list_struct = Vec::with_capacity(asps.num_ref_atlas_frame_lists_in_asps as usize);
        for _ in 0..asps.num_ref_atlas_frame_lists_in_asps {
            asps.ref_list_struct.push(RefListStruct::from_bitstream(
                bitstream,
                asps.long_term_ref_atlas_frames_flag,
                asps.log2_max_atlas_frame_order_cnt_lsb_minus_4 + 4,
            ));
        }
        asps.use_eight_orientations_flag = bitstream.read(1) != 0;
        asps.extended_projection_enabled_flag = bitstream.read(1) != 0;
        assert!(!asps.extended_projection_enabled_flag);
        if asps.extended_projection_enabled_flag {
            asps.max_number_projections_minus1 = bitstream.read_uvlc() as usize;
        }
        asps.normal_axis_limits_quantization_enabled_flag = bitstream.read(1) != 0;
        asps.normal_axis_max_delta_value_enabled_flag = bitstream.read(1) != 0;
        asps.patch_precedence_order_flag = bitstream.read(1) != 0;
        asps.log2_patch_packing_block_size = bitstream.read(3) as u8;
        asps.patch_size_quantizer_present_flag = bitstream.read(1) != 0;
        assert!(!asps.patch_size_quantizer_present_flag);
        asps.map_count_minus1 = bitstream.read(4) as u8;
        asps.pixel_deinterleaving_flag = bitstream.read(1) != 0;
        assert!(!asps.pixel_deinterleaving_flag);
        if asps.pixel_deinterleaving_flag {
            unimplemented!("pixel_deinterleaving_flag");
            // asps.pixel_deinterleaving_map_flag =
            //     Vec::with_capacity(asps.map_count_minus1 as usize + 1);
            // for _ in 0..asps.map_count_minus1 + 1 {
            //     asps.pixel_deinterleaving_map_flag
            //         .push(bitstream.read(1) != 0);
            // }
        }
        asps.raw_patch_enabled_flag = bitstream.read(1) != 0;
        asps.eom_patch_enabled_flag = bitstream.read(1) != 0;
        // (13Dec22)
        assert!(!asps.raw_patch_enabled_flag);
        assert!(!asps.eom_patch_enabled_flag);

        if asps.eom_patch_enabled_flag && asps.map_count_minus1 == 0 {
            asps.eom_fix_bit_count_minus1 = bitstream.read(4) as u8;
        }
        if asps.raw_patch_enabled_flag || asps.eom_patch_enabled_flag {
            asps.auxiliary_video_enabled_flag = bitstream.read(1) != 0;
        }
        // (13Dec22)
        assert!(!asps.auxiliary_video_enabled_flag);

        asps.plr_enabled_flag = bitstream.read(1) != 0;
        // (13Dec22)
        assert!(!asps.plr_enabled_flag);
        if asps.plr_enabled_flag {
            unimplemented!("PLR not implemented");
        }

        asps.vui_parameters_present_flag = bitstream.read(1) != 0;
        // (13Dec22)
        assert!(!asps.vui_parameters_present_flag);
        if asps.vui_parameters_present_flag {
            unimplemented!("VUI parameters are not implemented");
        }

        asps.extension_flag = bitstream.read(1) != 0;
        if asps.extension_flag {
            asps.vpcc_extension_flag = bitstream.read(1) != 0;
            asps.extension_7bits = bitstream.read(7) as u8;
        }

        if asps.vpcc_extension_flag {
            asps.vpcc_extension = AspsVpccExtension {
                remove_duplicate_point_enabled_flag: bitstream.read(1) != 0,
                surface_thickness_minus1: if asps.pixel_deinterleaving_flag || asps.plr_enabled_flag
                {
                    bitstream.read(7) as u8
                } else {
                    0
                },
            };
        }
        if asps.extension_7bits > 0 {
            unimplemented!("Extension 7 bits not implemented");
        }
        bitstream.byte_align();
        asps
    }
}

#[derive(Default, Debug, Clone)]
pub(crate) struct RefListStruct {
    pub(crate) num_ref_entries: u8,
    abs_delta_afoc_st: Vec<u8>,
    afoc_lsb_lt: Vec<u8>,
    st_ref_atlas_frame_flag: Vec<bool>,
    strpf_entry_sign_flag: Vec<bool>,
}

impl RefListStruct {
    fn from_bitstream(
        bitstream: &Bitstream,
        is_long_term_ref_atlas_frames: bool,
        log2_max_atlas_frame_order_cnt: u8,
    ) -> Self {
        let mut rls = RefListStruct::default();
        let num_entries = bitstream.read_uvlc() as usize;
        rls.num_ref_entries = num_entries as u8;

        // allocate space for the vectors
        rls.abs_delta_afoc_st = Vec::with_capacity(num_entries);
        rls.afoc_lsb_lt = Vec::with_capacity(num_entries);
        rls.st_ref_atlas_frame_flag = Vec::with_capacity(num_entries);
        rls.strpf_entry_sign_flag = Vec::with_capacity(num_entries);

        for _ in 0..rls.num_ref_entries {
            let is_st_ref_atlas_frame = if is_long_term_ref_atlas_frames {
                bitstream.read(1) != 0
            } else {
                true
            };
            rls.st_ref_atlas_frame_flag.push(is_st_ref_atlas_frame);

            if is_st_ref_atlas_frame {
                let abs_delta_afoc_st = bitstream.read_uvlc() as u8;
                rls.abs_delta_afoc_st.push(abs_delta_afoc_st);
                if abs_delta_afoc_st > 0 {
                    rls.strpf_entry_sign_flag.push(bitstream.read(1) != 0);
                } else {
                    rls.strpf_entry_sign_flag.push(true);
                }
            } else {
                rls.afoc_lsb_lt
                    .push(bitstream.read(log2_max_atlas_frame_order_cnt) as u8);
            }
        }
        rls
    }
}

#[derive(Default)]
pub(crate) struct AspsVpccExtension {
    pub(crate) remove_duplicate_point_enabled_flag: bool,
    pub(crate) surface_thickness_minus1: u8,
}

/// 8.3.6.2. Atlas Frame Parameter Set
#[derive(Default)]
pub(crate) struct AtlasFrameParameterSetRbsp {
    atlas_frame_parameter_set_id: u8,
    pub(crate) atlas_sequence_parameter_set_id: u8,
    pub(crate) atlas_frame_tile_information: AtlasFrameTileInformation,
    output_flag_present_flag: bool,
    pub(crate) num_ref_idx_default_active_minus1: u8,
    additional_lt_afoc_lsb_len: u8,
    lod_mode_enable_flag: bool,
    raw_3d_offset_bitcount_explicit_mode_flag: bool,
    extension_flag: bool,
    extension_8bits: u8,
    // afps_vpcc_extension: AfpsVpccExtension,
}

impl AtlasFrameParameterSetRbsp {
    fn from_bitstream(syntax: &Context, bitstream: &Bitstream) -> Self {
        let mut afps = AtlasFrameParameterSetRbsp::default();
        afps.atlas_frame_parameter_set_id = bitstream.read_uvlc() as u8;
        afps.atlas_sequence_parameter_set_id = bitstream.read_uvlc() as u8;
        afps.atlas_frame_tile_information = AtlasFrameTileInformation::from_bitstream(
            bitstream,
            syntax.get_atlas_sequence_parameter_set(afps.atlas_sequence_parameter_set_id as usize),
        );
        afps.output_flag_present_flag = bitstream.read(1) != 0;
        afps.num_ref_idx_default_active_minus1 = bitstream.read_uvlc() as u8;
        afps.additional_lt_afoc_lsb_len = bitstream.read_uvlc() as u8;
        afps.lod_mode_enable_flag = bitstream.read(1) != 0;
        afps.raw_3d_offset_bitcount_explicit_mode_flag = bitstream.read(1) != 0;
        afps.extension_flag = bitstream.read(1) != 0;
        if afps.extension_flag {
            afps.extension_8bits = bitstream.read(8) as u8;
        }
        if afps.extension_8bits > 0 {
            unimplemented!("afps extension 8 bits not implemented");
        }
        bitstream.byte_align();
        afps
    }
}

/// 8.3.6.2.2 Atlas Frame Tile Information Syntax
#[derive(Default)]
pub(crate) struct AtlasFrameTileInformation {
    /// If true, there is only 1 tile in each atlas frame (i.e. no partitioning)
    ///
    /// (12Dec22) Atlas Frame partition not implemented, i.e. this field is always true.
    pub single_tile_in_atlas_frame_flag: bool,
    uniform_partition_spacing_flag: bool,
    num_partition_columns_minus1: u32,
    num_partition_rows_minus1: u32,
    single_partition_per_tile_flag: u32,
    /// Since single_tile_in_atlas_frame is always true, this field is always 0
    pub num_tiles_in_atlas_frame_minus1: u32,
    /// (12Dec22) This field always returns false
    pub signalled_tile_id_flag: bool,
    signalled_tile_id_length_minus1: u8,
    /// if uniform_partition_spacing_flag is false, partition_column_width_minus1 will store the column width of columns 0..=num_partition_columns_minus1
    ///
    /// Used only in PCCDecoder::setTilePartitionSizeAfti
    partition_column_width_minus1: Vec<u32>,
    /// if uniform_partition_spacing_flag is false, partition_row_height_minus1 will store the row height of columns 0..=num_partition_rows_minus1
    ///
    /// Used only in PCCDecoder::setTilePartitionSizeAfti
    partition_row_height_minus1: Vec<u32>,
    top_left_partition_idx: Vec<u32>,
    bottom_right_partition_column_offset: Vec<u32>,
    bottom_right_partition_row_offset: Vec<u32>,
    // (13Dec22) tile_id is used only if signalled_tile_id_flag is true or in eom / raw patch data unit.
    // all of which is not supported for now.
    // Combined with num_tiles_in_atlas = 1, all references to tile_id should just return 0.
    // tile_id: Vec<u32>,
    auxiliary_video_tile_row_width_minus1: u32,
    auxiliary_video_tile_row_height: Vec<u32>,

    /// (12Dec22) Since single_tile_in_atlas_frame_flag is always true, these fields always have size 1.
    /// DIFF: changed to non-vec fields
    col_width: u16,
    row_height: u16,
    partition_pos: (usize, usize),
}

impl AtlasFrameTileInformation {
    fn from_bitstream(bitstream: &Bitstream, asps: &AtlasSequenceParameterSetRbsp) -> Self {
        let mut afti = AtlasFrameTileInformation {
            single_tile_in_atlas_frame_flag: bitstream.read(1) != 0,
            ..Default::default()
        };
        assert!(
            afti.single_tile_in_atlas_frame_flag,
            "(12Dec22) Atlas Frame partition not implemented"
        );
        if !afti.single_tile_in_atlas_frame_flag {
            unimplemented!("Atlas Frame partition not implemented");
        } else {
            afti.num_tiles_in_atlas_frame_minus1 = 0;
        }
        if asps.auxiliary_video_enabled_flag {
            afti.auxiliary_video_tile_row_width_minus1 = bitstream.read_uvlc();
            afti.auxiliary_video_tile_row_height
                .reserve_exact(afti.num_tiles_in_atlas_frame_minus1 as usize + 1);
            for _ in 0..(afti.num_tiles_in_atlas_frame_minus1 + 1) {
                afti.auxiliary_video_tile_row_height
                    .push(bitstream.read_uvlc());
            }
        }
        afti.signalled_tile_id_flag = bitstream.read(1) != 0;
        assert!(
            !afti.signalled_tile_id_flag,
            "(12Dec22) signalled_tile_id_flag always false"
        );
        // afti.tile_id
        //     .reserve_exact(afti.num_tiles_in_atlas_frame_minus1 as usize + 1);
        if afti.signalled_tile_id_flag {
            afti.signalled_tile_id_length_minus1 = bitstream.read_uvlc() as u8;
            // for _ in 0..=afti.num_tiles_in_atlas_frame_minus1 {
            //     afti.tile_id
            //         .push(bitstream.read(afti.signalled_tile_id_length_minus1 + 1));
            // }
            // afti.tile_id = v;
        } else {
            // see comments on afti.tile_id

            // for i in 0..=afti.num_tiles_in_atlas_frame_minus1 {
            //     afti.tile_id.push(i);
            // }
        }
        afti
    }

    #[inline]
    pub(crate) fn set_partition_width(&mut self, index: usize, width: u16) {
        assert!(index == 0);
        self.col_width = width;
    }

    #[inline]
    pub(crate) fn get_partition_width(&self, index: usize) -> u16 {
        assert!(index == 0);
        self.col_width
    }

    #[inline]
    pub(crate) fn set_partition_height(&mut self, index: usize, height: u16) {
        assert!(index == 0);
        self.row_height = height;
    }

    #[inline]
    pub(crate) fn get_partition_height(&self, index: usize) -> u16 {
        assert!(index == 0);
        self.row_height
    }
}

/// 8.3.6.4 Supplemental Enhancement Information (SEI) RBSP
/// Originally: PCCSEI
#[derive(Default)]
pub(crate) struct SeiRbsp {
    sei_prefix: Vec<Box<dyn SEI>>,
    sei_suffix: Vec<Box<dyn SEI>>,
}

#[derive(FromPrimitive, PartialEq, Clone, Copy, Debug)]
#[repr(u8)]
pub(crate) enum SeiPayloadType {
    #[default]
    BufferingPeriod = 0,
    AtlasFrameTiming,
    FillerPayload,
    UserDataRegisteredITUTT35,
    UserDataUnregistered,
    RecoveryPoint,
    NoReconstruction,
    TimeCode,
    SeiManifest,
    SeiPrefixIndication,
    ActiveSubBitstreams,
    ComponentCodecMapping,
    SceneObjectInformation,
    ObjectLabelInformation,
    PatchInformation,
    VolumetricRectangleInformation,
    AtlasObjectInformation,
    ViewportCameraParameters,
    ViewportPosition,
    DecodedAtlasInformationHash,
    AttributeTransformationParams = 64,
    OccupancySynthesis,
    GeometrySmoothing,
    AttributeSmoothing,
    ReservedSeiMessage,
}

impl SeiRbsp {
    fn from_bitstream(bitstream: &Bitstream, nal_unit_type: NalUnitType) -> Self {
        assert!(
            nal_unit_type == NalUnitType::PrefixESEI
                || nal_unit_type == NalUnitType::PrefixNSEI
                || nal_unit_type == NalUnitType::SuffixESEI
                || nal_unit_type == NalUnitType::SuffixNSEI
        );
        let mut payload_type = 0;
        loop {
            let byte = bitstream.read(8) as u8;
            payload_type += byte;
            if byte != 0xff {
                break;
            }
        }
        let payload_type = SeiPayloadType::from(payload_type);

        let mut payload_size = 0;
        loop {
            let byte = bitstream.read(8);
            payload_size += byte;
            if byte != 0xff {
                break;
            }
        }

        let mut sei_rbsp = Self {
            sei_prefix: Vec::new(),
            sei_suffix: Vec::new(),
        };

        if nal_unit_type == NalUnitType::PrefixESEI || nal_unit_type == NalUnitType::PrefixNSEI {
            match payload_type {
                SeiPayloadType::GeometrySmoothing => {
                    let sei = Box::new(SeiGeometrySmoothing::from_bitstream(bitstream));
                    sei_rbsp.sei_prefix.push(sei.clone());
                    sei
                }
                _ => unimplemented!("only geometry smoothing payload type is implemented"),
            };
        } else {
            unimplemented!("non-prefix-sei nal unit type not implemented");
        };

        bitstream.byte_align();
        // FIXME (3Dec22): somehow the tmc2 calls rbspTrailingBits after this which I can't find in the source code. super strange.
        // so added the following line which achieves the same effect as a hack.
        bitstream.read(8);
        sei_rbsp
    }

    pub(crate) fn is_sei_present(
        &self,
        nal_unit_type: NalUnitType,
        payload_type: SeiPayloadType,
    ) -> bool {
        if !nal_unit_type.is_prefix_sei() && !nal_unit_type.is_suffix_sei() {
            return false;
        }
        let sei = if nal_unit_type.is_prefix_sei() {
            &self.sei_prefix
        } else {
            &self.sei_suffix
        };
        sei.iter().any(|sei| sei.get_payload_type() == payload_type)
    }
}

/// Annex E: Supplemental Enhancement Information
/// F.2.1. General SEI message syntax  <=> 7.3.8 Supplemental enhancement Information message
trait SEI {
    fn get_payload_type(&self) -> SeiPayloadType;
}

/// H.20.2.19 Geometry smoothing SEI message syntax
#[derive(Default, Clone)]
struct SeiGeometrySmoothing {
    persistence_flag: bool,
    reset_flag: bool,
    instances_updated: u8,
    instance_index: Vec<u8>,
    instance_cancel_flag: Vec<bool>,
    method_type: Vec<u8>,
    filter_eom_points_flag: Vec<bool>,
    grid_size_minus_2: Vec<u8>,
    threshold: Vec<u8>,

    byte_str_data: Vec<u8>,
    payload_size: usize,
}

impl SeiGeometrySmoothing {
    fn from_bitstream(bitstream: &Bitstream) -> Self {
        let mut sei = SeiGeometrySmoothing::default();
        sei.persistence_flag = bitstream.read(1) != 0;
        sei.reset_flag = bitstream.read(1) != 0;
        sei.instances_updated = bitstream.read(8) as u8;

        sei.instance_index.resize(sei.instances_updated as usize, 0);
        sei.instance_cancel_flag
            .resize(sei.instances_updated as usize, false);
        sei.method_type.resize(sei.instances_updated as usize, 0);
        sei.filter_eom_points_flag
            .resize(sei.instances_updated as usize, false);
        sei.grid_size_minus_2
            .resize(sei.instances_updated as usize, 0);
        sei.threshold.resize(sei.instances_updated as usize, 0);

        for i in 0..sei.instances_updated as usize {
            sei.instance_index[i] = bitstream.read(8) as u8;
            let k = sei.instance_index[i] as usize;
            sei.instance_cancel_flag[k] = bitstream.read(1) != 0;
            if sei.instance_cancel_flag[k] {
                continue;
            }
            sei.method_type[k] = bitstream.read_uvlc() as u8;
            if sei.method_type[k] == 1 {
                sei.filter_eom_points_flag[k] = bitstream.read(1) != 0;
                sei.grid_size_minus_2[k] = bitstream.read(7) as u8;
                sei.threshold[k] = bitstream.read(8) as u8;
            }
        }

        sei
    }
}

impl SEI for SeiGeometrySmoothing {
    fn get_payload_type(&self) -> SeiPayloadType {
        SeiPayloadType::GeometrySmoothing
    }
}

/// 8.3.6.9  Atlas tile group layer (ATGL) Rbsp syntax
///
/// Stores the atlas data for each frame (and each tile in it)
pub(crate) struct AtlasTileLayerRbsp {
    pub(crate) header: AtlasTileHeader,
    /// assert!(self.tile_order == data_unit.tile_order)
    pub(crate) data_unit: AtlasTileDataUnit,
    pub(crate) atlas_frame_order_count_val: u32,
    pub(crate) atlas_frame_order_count_msb: u32,
    /// index in AtlasTileLayer
    // pub tile_order: usize,
    /// Frame Index (only used for encoding)
    pub(crate) enc_frame_index: usize,
    /// Tile Index (only used for encoding)
    pub(crate) enc_tile_index: usize,
    pub(crate) sei: Rc<Option<SeiRbsp>>,
}

impl AtlasTileLayerRbsp {
    fn from_bitstream(syntax: &Context, bitstream: &Bitstream, nal_unit_type: NalUnitType) -> Self {
        let header = AtlasTileHeader::from_bitstream(syntax, bitstream, nal_unit_type);
        let data_unit = AtlasTileDataUnit::from_bitstream(syntax, bitstream, &header);
        bitstream.byte_align();
        Self {
            header,
            data_unit,
            enc_frame_index: usize::MAX,
            enc_tile_index: usize::MAX,
            atlas_frame_order_count_msb: 0,
            atlas_frame_order_count_val: 0,
            sei: Rc::new(None),
            // FIXME: do we need to set tile order?
            // tile_order: 0,
        }
    }
}

/// 8.3.6.11  Atlas tile header syntax
#[derive(Default, Debug)]
pub(crate) struct AtlasTileHeader {
    no_output_of_prior_atlas_frames_flag: bool,
    pub(crate) frame_index: u8,
    pub(crate) atlas_frame_parameter_set_id: u8,
    atlas_adaptation_parameter_set_id: u8,
    pub(crate) id: u32,
    pub(crate) tile_type: TileType,
    atlas_output_flag: bool,
    pub(crate) atlas_frame_order_count_lsb: u32,
    pub(crate) ref_atlas_frame_list_sps_flag: bool,
    pub(crate) ref_atlas_frame_list_idx: u8,
    additional_afoc_lsb_present_flag: Vec<bool>,
    additional_afoc_lsb_val: Vec<u8>,
    pub(crate) pos_min_d_quantizer: u8,
    pos_delta_max_d_quantizer: u8,
    /// This is later copied to TileContext::log2_patch_quantizer
    pub(crate) patch_size_info_quantizer: (u8, u8),
    raw_3d_offset_axis_bitcount_minus1: u8,
    pub(crate) num_ref_idx_active_override_flag: bool,
    pub(crate) num_ref_idx_active_minus1: u8,
    pub(crate) ref_list_struct: RefListStruct,
    tile_nalu_type_info: u8,
}

impl AtlasTileHeader {
    fn from_bitstream(syntax: &Context, bitstream: &Bitstream, nal_unit_type: NalUnitType) -> Self {
        let mut ath = AtlasTileHeader::default();
        if nal_unit_type as u8 >= NalUnitType::BlaWLp as u8
            && nal_unit_type as u8 <= NalUnitType::Gcra as u8
        {
            ath.no_output_of_prior_atlas_frames_flag = bitstream.read(1) != 0;
        }

        if nal_unit_type == NalUnitType::TrailR {
            ath.tile_nalu_type_info = 1;
        }
        if nal_unit_type == NalUnitType::TrailN {
            ath.tile_nalu_type_info = 2;
        }

        ath.atlas_frame_parameter_set_id = bitstream.read_uvlc() as u8;
        ath.atlas_adaptation_parameter_set_id = bitstream.read_uvlc() as u8;
        let afps = syntax.get_atlas_frame_parameter_set(ath.atlas_frame_parameter_set_id as usize);
        let asps =
            syntax.get_atlas_sequence_parameter_set(afps.atlas_sequence_parameter_set_id as usize);
        let afti = &afps.atlas_frame_tile_information;

        if afti.signalled_tile_id_flag {
            ath.id = bitstream.read(afti.signalled_tile_id_length_minus1 + 1);
        } else if afti.num_tiles_in_atlas_frame_minus1 != 0 {
            let bit_count =
                fast_math::log2_raw(afti.num_tiles_in_atlas_frame_minus1 as f32 + 1.).ceil() as u8;
            ath.id = bitstream.read(bit_count);
        } else {
            ath.id = 0;
        }

        ath.tile_type = TileType::from(bitstream.read_uvlc() as u8);
        if afps.output_flag_present_flag {
            ath.atlas_output_flag = bitstream.read(1) != 0;
        } else {
            ath.atlas_output_flag = false;
        }

        ath.atlas_frame_order_count_lsb =
            bitstream.read(asps.log2_max_atlas_frame_order_cnt_lsb_minus_4 + 4);
        if asps.num_ref_atlas_frame_lists_in_asps > 0 {
            ath.ref_atlas_frame_list_sps_flag = bitstream.read(1) != 0;
        } else {
            ath.ref_atlas_frame_list_sps_flag = false;
        }

        ath.ref_atlas_frame_list_idx = 0;
        if !ath.ref_atlas_frame_list_sps_flag {
            ath.ref_list_struct = RefListStruct::from_bitstream(
                bitstream,
                asps.long_term_ref_atlas_frames_flag,
                asps.log2_max_atlas_frame_order_cnt_lsb_minus_4 + 4,
            );
        } else {
            ath.ref_list_struct =
                asps.ref_list_struct[ath.ref_atlas_frame_list_idx as usize].clone();
        }

        if asps.num_ref_atlas_frame_lists_in_asps > 1 {
            let bit_count =
                fast_math::log2_raw(asps.num_ref_atlas_frame_lists_in_asps as f32).ceil() as u8;
            ath.ref_atlas_frame_list_idx = bitstream.read(bit_count) as u8;
        }

        let rls_idx = ath.ref_atlas_frame_list_idx;
        let ref_list = if ath.ref_atlas_frame_list_sps_flag {
            &asps.ref_list_struct[rls_idx as usize]
        } else {
            &ath.ref_list_struct
        };

        let mut num_ltr_atlas_frame_entries = 0;
        for i in 0..ref_list.num_ref_entries {
            if !ref_list.st_ref_atlas_frame_flag[i as usize] {
                num_ltr_atlas_frame_entries += 1;
            }
        }

        for j in 0..num_ltr_atlas_frame_entries {
            ath.additional_afoc_lsb_present_flag
                .push(bitstream.read(1) != 0);
            if ath.additional_afoc_lsb_present_flag[j as usize] {
                ath.additional_afoc_lsb_val
                    .push(bitstream.read(afps.additional_lt_afoc_lsb_len) as u8);
            }
        }

        if ath.tile_type != TileType::Skip {
            if asps.normal_axis_limits_quantization_enabled_flag {
                ath.pos_min_d_quantizer = bitstream.read(5) as u8;
                if asps.normal_axis_limits_quantization_enabled_flag {
                    ath.pos_delta_max_d_quantizer = bitstream.read(5) as u8;
                }
            }
            if asps.patch_size_quantizer_present_flag {
                ath.patch_size_info_quantizer = (bitstream.read(3) as u8, bitstream.read(3) as u8);
            }
            if afps.raw_3d_offset_bitcount_explicit_mode_flag {
                let bit_count =
                    fast_math::log2_raw(asps.geometry_3d_bitdepth_minus1 as f32 + 1.).floor() as u8;
                ath.raw_3d_offset_axis_bitcount_minus1 = bitstream.read(bit_count) as u8;
            } else {
                // HMM possibly error.
                ath.raw_3d_offset_axis_bitcount_minus1 = std::cmp::max(
                    0,
                    asps.geometry_3d_bitdepth_minus1 - asps.geometry_2d_bitdepth_minus1,
                ) - 1;
            }
            if ath.tile_type == TileType::P && ref_list.num_ref_entries > 1 {
                ath.num_ref_idx_active_override_flag = bitstream.read(1) != 0;
                if ath.num_ref_idx_active_override_flag {
                    ath.num_ref_idx_active_minus1 = bitstream.read_uvlc() as u8;
                }
            }
        }
        bitstream.byte_align();
        ath
    }
}

#[derive(Default, Debug, PartialEq, FromPrimitive, Clone, Copy)]
#[repr(u8)]
pub(crate) enum TileType {
    /// Inter atlas tile
    #[default]
    P = 0,
    /// 1: Intra atlas tile
    I,
    /// 2: SKIP atlas tile
    Skip,
}

/// 8.3.7.1  General atlas tile data unit
pub(crate) struct AtlasTileDataUnit {
    /// assert!(self.tile_order == patch_information_data[i].tile_order)
    // tile_order: usize,
    pub(crate) patch_information_data: Vec<PatchInformationData>,
}

impl AtlasTileDataUnit {
    fn from_bitstream(syntax: &Context, bitstream: &Bitstream, ath: &AtlasTileHeader) -> Self {
        let mut atdu = AtlasTileDataUnit {
            // tile_order: 0,
            patch_information_data: Vec::new(),
        };

        if ath.tile_type == TileType::Skip {
            return atdu;
        }

        // let mut patch_index: usize = 0;
        // let prev_patch_size_u = 0;
        // let prev_patch_size_v = 0;
        // let pred_patch_index = 0;

        while let Some(pdu) = PatchInformationData::from_bitstream(syntax, bitstream, ath) {
            // FIXME
            // patch_index += 1
            atdu.patch_information_data.push(pdu);
        }

        // FIXME
        // prev_frame_index = atdu.tile_order;

        atdu
    }
}

#[derive(Debug, FromPrimitive)]
#[repr(u8)]
pub(crate) enum PatchModeITile {
    /// 0: Non-predicted patch mode
    #[default]
    Intra,
    /// 1: RAW point patch mode
    // Raw,
    /// 2: EOM (Enhanced Occupancy Mode) point patch mode
    /// Patch coding mode where a patch is associated with enhanced occupancy information
    // EOM,
    /// 14: Patch termination mode
    End = 14,
}

#[derive(Debug, FromPrimitive)]
#[repr(u8)]
pub(crate) enum PatchModePTile {
    #[default]
    /// 0: Patch skip mode
    Skip,
    /// 1: Patch merge mode
    Merge,
    /// 2: Inter predicted patch mode
    Inter,
    /// 3: Non-predicted patch mode
    Intra,
    /// 4: RAW point patch mode
    // Raw,
    /// 5: EOM point patch mode
    // EOM,
    /// 14: Patch termination mode
    End = 14,
}

/// 8.3.7.2  Patch information data syntax (pid)
#[derive(Debug)]
pub(crate) struct PatchInformationData {
    // tile_order: usize,
    // patch_index: usize,
    pub(crate) patch_mode: u8,
    pub(crate) patch_data_unit: PatchDataUnit,
}

#[derive(Debug)]
pub(crate) enum PatchDataUnit {
    Intra(IntraPatchDataUnit),
    Inter(InterPatchDataUnit),
    Merge(MergePatchDataUnit),
    Skip(SkipPatchDataUnit),
    // Raw(RawPatchDataUnit),
    // EOM(EOMPatchDataUnit),
}

impl PatchInformationData {
    /// Returns None if patchmode is PatchModePTile::End or PatchModeITile::End
    fn from_bitstream(
        syntax: &Context,
        bitstream: &Bitstream,
        ath: &AtlasTileHeader,
    ) -> Option<Self> {
        let patch_mode = bitstream.read_uvlc() as u8;
        match ath.tile_type {
            TileType::P => match PatchModePTile::from(patch_mode) {
                PatchModePTile::Merge => Some(PatchInformationData {
                    // tile_order: 0,
                    // patch_index: 0,
                    patch_mode,
                    patch_data_unit: PatchDataUnit::Merge(MergePatchDataUnit::from_bitstream(
                        syntax, bitstream, ath,
                    )),
                }),
                PatchModePTile::Inter => Some(PatchInformationData {
                    // tile_order: 0,
                    // patch_index: 0,
                    patch_mode,
                    patch_data_unit: PatchDataUnit::Inter(InterPatchDataUnit::from_bitstream(
                        syntax, bitstream, ath,
                    )),
                }),
                PatchModePTile::Intra => Some(PatchInformationData {
                    // tile_order: 0,
                    // patch_index: 0,
                    patch_mode,
                    patch_data_unit: PatchDataUnit::Intra(IntraPatchDataUnit::from_bitstream(
                        syntax, bitstream, ath,
                    )),
                }),
                PatchModePTile::Skip => Some(PatchInformationData {
                    // tile_order: 0,
                    // patch_index: 0,
                    patch_mode,
                    patch_data_unit: PatchDataUnit::Skip(SkipPatchDataUnit::default()),
                }),
                PatchModePTile::End => None,
            },
            TileType::I => match PatchModeITile::from(patch_mode) {
                PatchModeITile::Intra => Some(PatchInformationData {
                    // tile_order: 0,
                    // patch_index: 0,
                    patch_mode,
                    patch_data_unit: PatchDataUnit::Intra(IntraPatchDataUnit::from_bitstream(
                        syntax, bitstream, ath,
                    )),
                }),
                PatchModeITile::End => None,
            },
            TileType::Skip => unreachable!("TileType::Skip should not be used here"),
        }
    }
}

/// 8.3.7.3  Patch data unit syntax
/// Originally: PatchDataUnit. Used with patch_mode P_INTRA or I_INTRA
#[derive(Default, Debug)]
pub(crate) struct IntraPatchDataUnit {
    pub(crate) projection_id: u8,
    pub(crate) orientation_index: PatchOrientation,
    pub(crate) lod_enabled_flag: bool,
    // lod_scale_x_minus1: usize,
    // lod_scale_y_idc: u8,
    pub(crate) pos_2d: (usize, usize),         // (x,y)
    pub(crate) size_2d_minus1: (usize, usize), // (x,y)
    pub(crate) pos_3d_offset: (usize, usize),  // (u,v)
    pub(crate) pos_3d_offset_d: usize,
    pub(crate) pos_3d_range_d: usize,
    // point_local_reconstruction_data: PLRData,
    // patch_index: usize,
    // frame_index: usize,
    // tile_order: usize,
}

impl IntraPatchDataUnit {
    fn from_bitstream(syntax: &Context, bitstream: &Bitstream, ath: &AtlasTileHeader) -> Self {
        let mut pdu = IntraPatchDataUnit::default();
        let afps = syntax.get_atlas_frame_parameter_set(ath.atlas_frame_parameter_set_id as usize);
        let asps =
            syntax.get_atlas_sequence_parameter_set(afps.atlas_sequence_parameter_set_id as usize);
        let bitcount_uv = asps.geometry_3d_bitdepth_minus1 + 1;
        let bitcount_d = asps.geometry_3d_bitdepth_minus1 - ath.pos_min_d_quantizer + 1;

        pdu.pos_2d = (
            bitstream.read_uvlc() as usize,
            bitstream.read_uvlc() as usize,
        );
        pdu.size_2d_minus1 = (
            bitstream.read_uvlc() as usize,
            bitstream.read_uvlc() as usize,
        );
        pdu.pos_3d_offset = (
            bitstream.read(bitcount_uv) as usize,
            bitstream.read(bitcount_uv) as usize,
        );
        pdu.pos_3d_offset_d = bitstream.read(bitcount_d) as usize;

        if asps.normal_axis_max_delta_value_enabled_flag {
            let bitcount_for_max_depth = std::cmp::min(
                asps.geometry_2d_bitdepth_minus1,
                asps.geometry_3d_bitdepth_minus1,
            ) + 1
                - ath.pos_delta_max_d_quantizer;
            pdu.pos_3d_range_d = bitstream.read(bitcount_for_max_depth) as usize;
        }

        pdu.projection_id = bitstream
            .read(fast_math::log2_raw(asps.max_number_projections_minus1 as f32 + 1.).ceil() as u8)
            as u8;
        assert!(pdu.projection_id <= 5);
        pdu.orientation_index =
            PatchOrientation::from(bitstream.read(if asps.use_eight_orientations_flag {
                3
            } else {
                1
            }) as u8);

        if afps.lod_mode_enable_flag {
            unimplemented!("lod_mode_enable_flag")
        }
        if asps.plr_enabled_flag {
            unimplemented!("plr_enabled_flag")
        }
        pdu
    }
}

#[derive(Default, Debug)]
pub(crate) struct InterPatchDataUnit {
    ref_index: usize,
    ref_patch_index: usize,

    pos_2d: (i32, i32),        // (x, y)
    delta_2d_size: (i32, i32), // (x, y)
    pos_3d_offset: (i32, i32), // (u, v)
    pos_3d_offset_d: i32,
    pos_3d_range_d: i32,
    // point_local_reconstruction_data: PLRData,
    // patch_index: usize,
    // frame_index: usize,
    // tile_order: usize,
}

impl InterPatchDataUnit {
    fn from_bitstream(syntax: &Context, bitstream: &Bitstream, ath: &AtlasTileHeader) -> Self {
        let afps = syntax.get_atlas_frame_parameter_set(ath.atlas_frame_parameter_set_id as usize);
        let asps =
            syntax.get_atlas_sequence_parameter_set(afps.atlas_sequence_parameter_set_id as usize);
        let num_ref_idx_active = syntax.get_num_ref_idx_active(ath);

        let mut pdu = InterPatchDataUnit {
            ref_index: if num_ref_idx_active > 1 {
                bitstream.read_uvlc() as usize
            } else {
                0
            },
            ref_patch_index: bitstream.read_svlc() as usize,
            pos_2d: (bitstream.read_svlc(), bitstream.read_svlc()),
            delta_2d_size: (bitstream.read_svlc(), bitstream.read_svlc()),
            pos_3d_offset: (bitstream.read_svlc(), bitstream.read_svlc()),
            pos_3d_offset_d: bitstream.read_svlc(),

            ..Default::default()
        };

        if asps.normal_axis_max_delta_value_enabled_flag {
            unimplemented!("normal_axis_max_delta_value_enabled_flag")
        }

        if asps.plr_enabled_flag {
            unimplemented!("plr_enabled_flag")
        }
        pdu
    }
}

#[derive(Default, Debug)]
pub(crate) struct MergePatchDataUnit {
    override_2d_params_flag: bool,
    override_3d_params_flag: bool,
    override_plr_flag: bool,
    ref_index: usize,

    pos_2d_x: i32,
    pos_2d_y: i32,
    delta_2d_size_x: i32,
    delta_2d_size_y: i32,
    pos_3d_offset_u: i32,
    pos_3d_offset_v: i32,
    pos_3d_offset_d: i32,
    pos_3d_range_d: i32,
    // point_local_reconstruction_data: PLRData,
    // patch_index: usize,
    // frame_index: usize,
    // tile_order: usize,
}

impl MergePatchDataUnit {
    fn from_bitstream(syntax: &Context, bitstream: &Bitstream, ath: &AtlasTileHeader) -> Self {
        let mut pdu = MergePatchDataUnit::default();
        let afps = syntax.get_atlas_frame_parameter_set(ath.atlas_frame_parameter_set_id as usize);
        let asps =
            syntax.get_atlas_sequence_parameter_set(afps.atlas_sequence_parameter_set_id as usize);
        let num_ref_idx_active = syntax.get_num_ref_idx_active(ath);

        if num_ref_idx_active > 1 {
            pdu.ref_index = bitstream.read_uvlc() as usize;
        }
        pdu.override_2d_params_flag = bitstream.read(1) != 0;

        if pdu.override_2d_params_flag {
            pdu.pos_2d_x = bitstream.read_svlc();
            pdu.pos_2d_y = bitstream.read_svlc();
            pdu.delta_2d_size_x = bitstream.read_svlc();
            pdu.delta_2d_size_y = bitstream.read_svlc();

            if asps.plr_enabled_flag {
                unimplemented!("plr_enabled_flag")
            }
        } else {
            pdu.override_3d_params_flag = bitstream.read(1) != 0;
            pdu.pos_3d_offset_u = bitstream.read_svlc();
            pdu.pos_3d_offset_v = bitstream.read_svlc();
            pdu.pos_3d_offset_d = bitstream.read_svlc();
            if asps.normal_axis_max_delta_value_enabled_flag {
                // pdu.pos_3d_range_d = bitstream.read_svlc();
                unimplemented!("normal_axis_max_delta_value_enabled_flag")
            }

            if asps.plr_enabled_flag {
                unimplemented!("plr_enabled_flag")
            }
        }

        if asps.plr_enabled_flag {
            unimplemented!("plr_enabled_flag")
        }
        pdu
    }
}

#[derive(Default, Debug)]
pub(crate) struct SkipPatchDataUnit {}

// #[derive(Default, Debug)]
// struct RawPatchDataUnit {
//     patch_in_auxiliary_video_flag: bool,
//     pos_2d_x: usize,
//     pos_2d_y: usize,
//     size_2d_x_minus1: usize,
//     size_2d_y_minus1: usize,
//     pos_3d_offset_u: usize,
//     pos_3d_offset_v: usize,
//     pos_3d_offset_d: usize,
//     pos_3d_range_d: usize,
//     raw_points_minus1: u32,
//     // patch_index: usize,
//     // frame_index: usize,
//     // tile_order: usize,
// }

// Enhanced Occupancy Mode (EOM)
//
// #[derive(Default, Debug)]
// struct EOMPatchDataUnit {
//     patch_in_auxiliary_video_flag: bool,

//     pos_2d_x: usize,
//     pos_2d_y: usize,
//     size_2d_x_minus1: usize,
//     size_2d_y_minus1: usize,

//     patch_count_minus1: usize,
//     associated_patches_idx: Vec<usize>,
//     points: Vec<usize>,
//     // patch_index: usize,
//     // frame_index: usize,
//     // tile_order: usize,
// }

// pub struct Reader {
//     prev_patch_size_u: i32,
//     prev_patch_size_v: i32,
//     pred_patch_idx: i32,
//     prev_frame_idx: i32,
// }

// impl Reader {}
