//! Point-in-polygon (PIP) tests.
//!
//! Two classic algorithms are implemented on rings (closed polylines):
//!
//! * **Ray casting / crossing number** ([`ring_contains_crossing`]): shoot a
//!   horizontal ray from the query point towards `+x` and count how many edges
//!   it crosses. Odd means inside. This is the **even-odd** fill rule.
//! * **Winding number** ([`ring_winding_number`], Dan Sunday's formulation):
//!   count how many times the ring winds around the point, `+1` for each
//!   upward edge crossing the ray with the point on its left, `-1` for each
//!   downward one with the point on its right. Non-zero means inside. This is
//!   the **non-zero** fill rule.
//!
//! For simple (non self-intersecting) rings both rules agree. They differ for
//! self-intersecting rings, e.g. the centre of a pentagram is "outside" under
//! even-odd (crossed twice) but "inside" under non-zero (wound twice).
//!
//! ## Rings
//!
//! A ring is a `&[Point]`. It may be given open (`[a, b, c]`) or closed
//! (`[a, b, c, a]`); the closing edge from the last point back to the first is
//! always implied, and a repeated closing point just adds a zero-length edge,
//! which both algorithms ignore. Orientation (CW or CCW) does not matter.
//!
//! ## Boundary behaviour (the contract)
//!
//! The raw ring functions ([`ring_contains_crossing`], [`ring_winding_number`])
//! use the standard *half-open* rule: an edge counts if one endpoint is
//! strictly above the point and the other is on or below it. That makes every
//! point classify consistently (a vertex is never counted twice), but for
//! points exactly **on** the boundary the answer is arbitrary-but-deterministic
//! (typically: left/bottom edges are "in", right/top edges are "out").
//!
//! [`Polygon::locate`] resolves this explicitly: it first checks whether the
//! point lies on any ring (exterior or hole) and, if so, returns
//! [`Location::Boundary`]. The on-segment test is **exact in floating point**
//! (`orient == 0.0` plus a bounding-box check), with no epsilon. A point
//! computed by interpolation along a slanted edge may therefore land a hair
//! off the edge and be reported as `Inside` or `Outside`; axis-aligned edges
//! and vertices are detected exactly. [`Polygon::contains`] is strict: it
//! returns `true` only for [`Location::Inside`].
//!
//! ## Holes
//!
//! A [`Polygon`] has one exterior ring and zero or more interior rings
//! (holes). A point is inside the polygon when it is inside the exterior and
//! not inside any hole. Holes are assumed to lie within the exterior and not
//! to overlap each other (the usual OGC "valid polygon" assumptions); this is
//! not checked.

use crate::types::{orient, Point, Rect};

/// Result of locating a point relative to a polygon.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Location {
    /// Strictly inside the polygon's area (not on any ring).
    Inside,
    /// Strictly outside (including inside a hole, off its boundary).
    Outside,
    /// Exactly on the exterior ring or on a hole's ring.
    Boundary,
}

/// Which algorithm (and therefore which fill rule) to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum PipAlgorithm {
    /// Crossing number, even-odd fill rule.
    RayCasting,
    /// Winding number, non-zero fill rule.
    #[default]
    WindingNumber,
}

/// Iterates over the edges `(p[i], p[i+1])` of a ring, including the closing
/// edge `(p[n-1], p[0])`.
fn edges(ring: &[Point]) -> impl Iterator<Item = (Point, Point)> + '_ {
    ring.iter()
        .copied()
        .zip(ring.iter().copied().cycle().skip(1))
}

/// `true` if `p` lies exactly on the closed segment `[a, b]`.
pub fn point_on_segment(p: Point, a: Point, b: Point) -> bool {
    orient(a, b, p) == 0.0
        && p.x >= a.x.min(b.x)
        && p.x <= a.x.max(b.x)
        && p.y >= a.y.min(b.y)
        && p.y <= a.y.max(b.y)
}

/// `true` if `p` lies exactly on one of the ring's edges.
pub fn point_on_ring(ring: &[Point], p: Point) -> bool {
    edges(ring).any(|(a, b)| point_on_segment(p, a, b))
}

/// Ray casting (crossing number, even-odd rule) for a single ring.
///
/// See the module docs for what happens on the boundary.
pub fn ring_contains_crossing(ring: &[Point], p: Point) -> bool {
    let mut inside = false;
    for (a, b) in edges(ring) {
        // Half-open rule: the edge straddles the horizontal line y = p.y.
        // This also skips horizontal edges (a.y == b.y), avoiding a division by zero.
        if (a.y > p.y) != (b.y > p.y) {
            // x coordinate where the edge crosses y = p.y.
            let x_cross = a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y);
            if p.x < x_cross {
                inside = !inside;
            }
        }
    }
    inside
}

/// Winding number of a single ring around `p` (Sunday's algorithm).
///
/// Positive for counter-clockwise rings, negative for clockwise ones, zero
/// when `p` is outside. Only `orient` is used: no division, no trigonometry.
pub fn ring_winding_number(ring: &[Point], p: Point) -> i32 {
    let mut wn = 0;
    for (a, b) in edges(ring) {
        if a.y <= p.y {
            // Upward crossing with p strictly left of a->b.
            if b.y > p.y && orient(a, b, p) > 0.0 {
                wn += 1;
            }
        } else if b.y <= p.y && orient(a, b, p) < 0.0 {
            // Downward crossing with p strictly right of a->b.
            wn -= 1;
        }
    }
    wn
}

/// Tests a single ring with the chosen algorithm (boundary is arbitrary, see
/// the module docs).
pub fn ring_contains(ring: &[Point], p: Point, algorithm: PipAlgorithm) -> bool {
    match algorithm {
        PipAlgorithm::RayCasting => ring_contains_crossing(ring, p),
        PipAlgorithm::WindingNumber => ring_winding_number(ring, p) != 0,
    }
}

/// A polygon with an exterior ring and optional holes.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Polygon {
    pub exterior: Vec<Point>,
    pub interiors: Vec<Vec<Point>>,
}

impl Polygon {
    pub fn new(exterior: Vec<Point>, interiors: Vec<Vec<Point>>) -> Self {
        Self {
            exterior,
            interiors,
        }
    }

    /// Polygon without holes.
    pub fn without_holes(exterior: Vec<Point>) -> Self {
        Self::new(exterior, Vec::new())
    }

    /// Bounding box of the exterior ring.
    pub fn bbox(&self) -> Rect {
        Rect::from_points(&self.exterior)
    }

    /// Iterates over all rings, exterior first.
    pub fn rings(&self) -> impl Iterator<Item = &[Point]> {
        std::iter::once(self.exterior.as_slice()).chain(self.interiors.iter().map(Vec::as_slice))
    }

    /// Locates `p` using the winding-number algorithm.
    pub fn locate(&self, p: Point) -> Location {
        self.locate_with(p, PipAlgorithm::WindingNumber)
    }

    /// Locates `p` with an explicit algorithm. Boundary points are reported
    /// as [`Location::Boundary`] regardless of the algorithm.
    pub fn locate_with(&self, p: Point, algorithm: PipAlgorithm) -> Location {
        if self.exterior.is_empty() || !self.bbox().contains_point(p) {
            return Location::Outside;
        }
        if self.rings().any(|ring| point_on_ring(ring, p)) {
            return Location::Boundary;
        }
        if !ring_contains(&self.exterior, p, algorithm) {
            return Location::Outside;
        }
        if self
            .interiors
            .iter()
            .any(|hole| ring_contains(hole, p, algorithm))
        {
            return Location::Outside;
        }
        Location::Inside
    }

    /// `true` only when `p` is strictly inside (boundary excluded).
    pub fn contains(&self, p: Point) -> bool {
        self.locate(p) == Location::Inside
    }

    /// `true` when `p` is inside or on the boundary.
    pub fn intersects_point(&self, p: Point) -> bool {
        self.locate(p) != Location::Outside
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use PipAlgorithm::*;

    fn pts(v: &[(f64, f64)]) -> Vec<Point> {
        v.iter().map(|&(x, y)| Point::new(x, y)).collect()
    }

    fn square(min: f64, max: f64) -> Vec<Point> {
        pts(&[(min, min), (max, min), (max, max), (min, max)])
    }

    #[test]
    fn simple_square_both_algorithms() {
        let poly = Polygon::without_holes(square(0.0, 10.0));
        for alg in [RayCasting, WindingNumber] {
            assert_eq!(poly.locate_with(Point::new(5.0, 5.0), alg), Location::Inside);
            assert_eq!(poly.locate_with(Point::new(15.0, 5.0), alg), Location::Outside);
            assert_eq!(poly.locate_with(Point::new(-1.0, 5.0), alg), Location::Outside);
            assert_eq!(poly.locate_with(Point::new(5.0, 11.0), alg), Location::Outside);
        }
    }

    #[test]
    fn boundary_is_detected() {
        let poly = Polygon::without_holes(square(0.0, 10.0));
        for alg in [RayCasting, WindingNumber] {
            for p in [(0.0, 0.0), (10.0, 5.0), (5.0, 10.0), (0.0, 7.5), (10.0, 10.0)] {
                assert_eq!(
                    poly.locate_with(Point::new(p.0, p.1), alg),
                    Location::Boundary,
                    "{p:?} {alg:?}"
                );
            }
        }
        assert!(!poly.contains(Point::new(0.0, 5.0)));
        assert!(poly.intersects_point(Point::new(0.0, 5.0)));
    }

    #[test]
    fn raw_ring_half_open_rule_is_consistent() {
        // The raw ring tests are deterministic on the boundary: the left edge
        // counts as inside and the right edge as outside for this CCW square.
        let ring = square(0.0, 10.0);
        assert!(ring_contains_crossing(&ring, Point::new(0.0, 5.0)));
        assert!(!ring_contains_crossing(&ring, Point::new(10.0, 5.0)));
        assert_eq!(ring_winding_number(&ring, Point::new(0.0, 5.0)), 1);
        assert_eq!(ring_winding_number(&ring, Point::new(10.0, 5.0)), 0);
    }

    #[test]
    fn slanted_edge_boundary() {
        let tri = Polygon::without_holes(pts(&[(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)]));
        assert_eq!(tri.locate(Point::new(2.0, 2.0)), Location::Boundary);
        assert_eq!(tri.locate(Point::new(1.0, 1.0)), Location::Inside);
        assert_eq!(tri.locate(Point::new(2.5, 2.5)), Location::Outside);
    }

    #[test]
    fn orientation_does_not_matter() {
        let mut cw = square(0.0, 10.0);
        cw.reverse();
        let p = Point::new(3.0, 3.0);
        assert_eq!(ring_winding_number(&cw, p), -1);
        assert_eq!(ring_winding_number(&square(0.0, 10.0), p), 1);
        assert!(Polygon::without_holes(cw).contains(p));
    }

    #[test]
    fn closed_ring_repeated_point_is_fine() {
        let mut ring = square(0.0, 10.0);
        ring.push(ring[0]);
        let poly = Polygon::without_holes(ring);
        assert!(poly.contains(Point::new(1.0, 1.0)));
        assert!(!poly.contains(Point::new(11.0, 1.0)));
    }

    #[test]
    fn polygon_with_hole() {
        let poly = Polygon::new(square(0.0, 10.0), vec![square(4.0, 6.0)]);
        for alg in [RayCasting, WindingNumber] {
            assert_eq!(poly.locate_with(Point::new(2.0, 2.0), alg), Location::Inside);
            assert_eq!(poly.locate_with(Point::new(5.0, 5.0), alg), Location::Outside);
            assert_eq!(poly.locate_with(Point::new(4.0, 5.0), alg), Location::Boundary);
            assert_eq!(poly.locate_with(Point::new(6.0, 6.0), alg), Location::Boundary);
        }
    }

    #[test]
    fn concave_polygon() {
        // A "U" shape: the notch (5, 8) is outside.
        let u = Polygon::without_holes(pts(&[
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (7.0, 10.0),
            (7.0, 3.0),
            (3.0, 3.0),
            (3.0, 10.0),
            (0.0, 10.0),
        ]));
        for alg in [RayCasting, WindingNumber] {
            assert_eq!(u.locate_with(Point::new(5.0, 8.0), alg), Location::Outside);
            assert_eq!(u.locate_with(Point::new(5.0, 1.0), alg), Location::Inside);
            assert_eq!(u.locate_with(Point::new(1.0, 9.0), alg), Location::Inside);
            // The ray from here runs along the horizontal edge (3,3)-(7,3).
            assert_eq!(u.locate_with(Point::new(1.0, 3.0), alg), Location::Inside);
        }
    }

    #[test]
    fn ray_through_vertex_is_counted_once() {
        // Diamond: a ray from the centre-left passes exactly through vertex (10, 5).
        let diamond = pts(&[(5.0, 0.0), (10.0, 5.0), (5.0, 10.0), (0.0, 5.0)]);
        let p = Point::new(2.0, 5.0);
        assert!(ring_contains_crossing(&diamond, p));
        assert_eq!(ring_winding_number(&diamond, p), 1);
        let outside = Point::new(-2.0, 5.0);
        assert!(Polygon::without_holes(diamond).locate(outside) == Location::Outside);
    }

    #[test]
    fn pentagram_even_odd_vs_non_zero() {
        // Self-intersecting 5-point star drawn in one stroke.
        let star: Vec<Point> = (0..5)
            .map(|i| {
                let a = std::f64::consts::FRAC_PI_2 + (i as f64) * 4.0 * std::f64::consts::PI / 5.0;
                Point::new(a.cos() * 10.0, a.sin() * 10.0)
            })
            .collect();
        let poly = Polygon::without_holes(star);
        let centre = Point::new(0.0, 0.0);
        assert_eq!(poly.locate_with(centre, RayCasting), Location::Outside);
        assert_eq!(poly.locate_with(centre, WindingNumber), Location::Inside);
        assert_eq!(ring_winding_number(&poly.exterior, centre).abs(), 2);
    }

    #[test]
    fn degenerate_inputs() {
        let empty = Polygon::default();
        assert_eq!(empty.locate(Point::new(0.0, 0.0)), Location::Outside);
        let line = Polygon::without_holes(pts(&[(0.0, 0.0), (1.0, 1.0)]));
        assert_eq!(line.locate(Point::new(0.5, 0.5)), Location::Boundary);
        assert_eq!(line.locate(Point::new(0.5, 0.4)), Location::Outside);
    }
}
