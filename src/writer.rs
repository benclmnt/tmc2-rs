use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use crate::codec::PointSet3;

// http://gamma.cs.unc.edu/POWERPLANT/papers/ply.pdf
pub enum Format {
    Ascii,
    // BinaryLittleEndian,
    // BinaryBigEndian,
}

pub(crate) struct PlyWriter {
    pc: PointSet3,
    format: Format,
}

impl PlyWriter {
    pub fn new(pc: PointSet3, format: Format) -> Self {
        Self { pc, format }
    }

    pub fn write(&self, path: &Path) {
        let output = File::create(path).unwrap();
        let mut w = BufWriter::new(output);
        self.write_header(&mut w);
        self.write_body(&mut w);
    }

    fn write_header<P>(&self, file: &mut P)
    where
        P: std::io::Write,
    {
        writeln!(file, "ply").unwrap();
        match self.format {
            Format::Ascii => {
                writeln!(file, "format ascii 1.0").unwrap();
            } // Format::BinaryLittleEndian => {
              //     writeln!(file, "format binary_little_endian 1.0").unwrap();
              // }
              // Format::BinaryBigEndian => {
              //     writeln!(file, "format binary_big_endian 1.0").unwrap();
              // }
        }
        writeln!(file, "element vertex {}", self.pc.point_count()).unwrap();
        writeln!(file, "property uint x").unwrap();
        writeln!(file, "property uint y").unwrap();
        writeln!(file, "property uint z").unwrap();
        if self.pc.with_colors {
            writeln!(file, "property uchar red").unwrap();
            writeln!(file, "property uchar green").unwrap();
            writeln!(file, "property uchar blue").unwrap();
        }
        writeln!(file, "element face 0").unwrap();
        writeln!(file, "property list uint8 int32 vertex_index").unwrap();
        writeln!(file, "end_header").unwrap();
    }

    fn write_body<P>(&self, file: &mut P)
    where
        P: std::io::Write,
    {
        for i in 0..self.pc.point_count() {
            let pos = &self.pc.positions[i];
            write!(file, "{} {} {}", pos.x, pos.y, pos.z).unwrap();
            if self.pc.with_colors {
                let color = &self.pc.colors[i];
                write!(file, " {} {} {}", color.x, color.y, color.z).unwrap();
            }
            writeln!(file).unwrap();
        }
    }
}
