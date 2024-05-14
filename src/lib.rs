//! # geo-algos
//!
//! Small, dependency-free implementations of classic geospatial algorithms,
//! written from scratch for learning purposes.
//!
//! | Module        | What it contains                                                    |
//! |---------------|---------------------------------------------------------------------|
//! | [`types`]     | [`Point`] (planar x/y), [`LonLat`] (degrees), [`Rect`] (AABB)       |
//! | [`pip`]       | Point-in-polygon: ray casting (even-odd) and winding number         |
//! | [`distance`]  | Haversine (sphere) and Vincenty inverse (WGS84 ellipsoid)           |
//! | [`simplify`]  | Douglas-Peucker polyline simplification                             |
//! | [`rtree`]     | R-tree: insert (linear/quadratic split), STR bulk load, range, kNN  |
//!
//! Everything works on `f64`. There is no exact/robust arithmetic: see the
//! README for the consequences.

pub mod distance;
pub mod pip;
pub mod rtree;
pub mod simplify;
pub mod types;

pub use distance::{haversine, vincenty_inverse, Geodesic, VincentyError};
pub use pip::{Location, PipAlgorithm, Polygon};
pub use rtree::{Entry, RTree, SplitStrategy};
pub use simplify::{douglas_peucker, douglas_peucker_indices};
pub use types::{LonLat, Point, Rect};
