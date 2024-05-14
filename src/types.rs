//! Basic geometric types shared by every module.

/// A point in a planar (Cartesian) coordinate system.
///
/// Planar algorithms (point-in-polygon, Douglas-Peucker, R-tree) use this
/// type. If you feed them longitude/latitude, they treat the degrees as flat
/// x/y coordinates (an equirectangular view), which is fine for indexing and
/// small areas but is not a geodesic computation.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Squared Euclidean distance (cheaper: no square root).
    pub fn distance_sq(self, other: Point) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    /// Euclidean distance.
    pub fn distance(self, other: Point) -> f64 {
        self.distance_sq(other).sqrt()
    }
}

impl From<(f64, f64)> for Point {
    fn from((x, y): (f64, f64)) -> Self {
        Point::new(x, y)
    }
}

impl From<LonLat> for Point {
    /// Maps longitude to `x` and latitude to `y` (degrees, unprojected).
    fn from(ll: LonLat) -> Self {
        Point::new(ll.lon, ll.lat)
    }
}

/// Twice the signed area of the triangle `(a, b, c)`.
///
/// * `> 0`: `c` is to the left of the directed line `a -> b` (counter-clockwise turn)
/// * `< 0`: `c` is to the right (clockwise turn)
/// * `= 0`: the three points are collinear (in floating point)
pub fn orient(a: Point, b: Point, c: Point) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)
}

/// A geographic position in **degrees**: longitude first, like GeoJSON.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct LonLat {
    /// Longitude in degrees, east positive, expected in `[-180, 180]`.
    pub lon: f64,
    /// Latitude in degrees, north positive, expected in `[-90, 90]`.
    pub lat: f64,
}

impl LonLat {
    pub const fn new(lon: f64, lat: f64) -> Self {
        Self { lon, lat }
    }

    /// Converts degrees/minutes/seconds to decimal degrees. The sign is taken
    /// from `deg` (use `-0.0` for values between -1 and 0 degrees).
    pub fn dms_to_deg(deg: f64, min: f64, sec: f64) -> f64 {
        let magnitude = deg.abs() + min / 60.0 + sec / 3600.0;
        if deg.is_sign_negative() {
            -magnitude
        } else {
            magnitude
        }
    }

    /// `true` when both coordinates are finite and in the usual ranges.
    pub fn is_valid(&self) -> bool {
        self.lon.is_finite()
            && self.lat.is_finite()
            && (-180.0..=180.0).contains(&self.lon)
            && (-90.0..=90.0).contains(&self.lat)
    }
}

/// An axis-aligned rectangle (bounding box). Boundaries are **closed**:
/// a point on an edge is contained, and rectangles that only touch intersect.
///
/// [`Rect::empty`] is the identity for [`Rect::union`]: its `min` is `+inf`
/// and its `max` is `-inf`, so it contains and intersects nothing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rect {
    pub min: Point,
    pub max: Point,
}

impl Rect {
    /// Rectangle spanned by two corner points (in any order).
    pub fn new(a: Point, b: Point) -> Self {
        Self {
            min: Point::new(a.x.min(b.x), a.y.min(b.y)),
            max: Point::new(a.x.max(b.x), a.y.max(b.y)),
        }
    }

    /// Degenerate rectangle containing a single point.
    pub fn from_point(p: Point) -> Self {
        Self { min: p, max: p }
    }

    /// Smallest rectangle containing all the points (empty if none).
    pub fn from_points<'a>(points: impl IntoIterator<Item = &'a Point>) -> Self {
        points
            .into_iter()
            .fold(Rect::empty(), |r, p| r.union(&Rect::from_point(*p)))
    }

    /// The empty rectangle (identity element for [`Rect::union`]).
    pub fn empty() -> Self {
        Self {
            min: Point::new(f64::INFINITY, f64::INFINITY),
            max: Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        }
    }

    pub fn is_empty(&self) -> bool {
        !(self.min.x <= self.max.x && self.min.y <= self.max.y)
    }

    pub fn width(&self) -> f64 {
        if self.is_empty() {
            0.0
        } else {
            self.max.x - self.min.x
        }
    }

    pub fn height(&self) -> f64 {
        if self.is_empty() {
            0.0
        } else {
            self.max.y - self.min.y
        }
    }

    pub fn area(&self) -> f64 {
        self.width() * self.height()
    }

    /// Half-perimeter (called "margin" in the R*-tree paper).
    pub fn margin(&self) -> f64 {
        self.width() + self.height()
    }

    pub fn center(&self) -> Point {
        Point::new(
            (self.min.x + self.max.x) * 0.5,
            (self.min.y + self.max.y) * 0.5,
        )
    }

    /// Smallest rectangle containing both.
    pub fn union(&self, other: &Rect) -> Rect {
        Rect {
            min: Point::new(self.min.x.min(other.min.x), self.min.y.min(other.min.y)),
            max: Point::new(self.max.x.max(other.max.x), self.max.y.max(other.max.y)),
        }
    }

    /// How much the area grows if `other` is added: `area(self ∪ other) - area(self)`.
    pub fn enlargement(&self, other: &Rect) -> f64 {
        self.union(other).area() - self.area()
    }

    /// Closed intersection test (touching rectangles intersect).
    pub fn intersects(&self, other: &Rect) -> bool {
        self.min.x <= other.max.x
            && other.min.x <= self.max.x
            && self.min.y <= other.max.y
            && other.min.y <= self.max.y
    }

    pub fn contains_point(&self, p: Point) -> bool {
        self.min.x <= p.x && p.x <= self.max.x && self.min.y <= p.y && p.y <= self.max.y
    }

    pub fn contains_rect(&self, other: &Rect) -> bool {
        other.is_empty()
            || (self.min.x <= other.min.x
                && other.max.x <= self.max.x
                && self.min.y <= other.min.y
                && other.max.y <= self.max.y)
    }

    /// Squared distance from `p` to the closest point of the rectangle
    /// (0 if `p` is inside). This is the "MINDIST" used by kNN search: a
    /// lower bound of the distance from `p` to anything stored inside.
    pub fn min_distance_sq(&self, p: Point) -> f64 {
        if self.is_empty() {
            return f64::INFINITY;
        }
        let dx = (self.min.x - p.x).max(0.0).max(p.x - self.max.x);
        let dy = (self.min.y - p.y).max(0.0).max(p.y - self.max.y);
        dx * dx + dy * dy
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orient_signs() {
        let a = Point::new(0.0, 0.0);
        let b = Point::new(1.0, 0.0);
        assert!(orient(a, b, Point::new(0.5, 1.0)) > 0.0);
        assert!(orient(a, b, Point::new(0.5, -1.0)) < 0.0);
        assert_eq!(orient(a, b, Point::new(2.0, 0.0)), 0.0);
    }

    #[test]
    fn rect_basics() {
        let r = Rect::new(Point::new(2.0, 3.0), Point::new(0.0, 1.0));
        assert_eq!(r.min, Point::new(0.0, 1.0));
        assert_eq!(r.max, Point::new(2.0, 3.0));
        assert_eq!(r.area(), 4.0);
        assert_eq!(r.margin(), 4.0);
        assert_eq!(r.center(), Point::new(1.0, 2.0));
        assert!(r.contains_point(Point::new(2.0, 3.0))); // closed boundary
        assert!(!r.contains_point(Point::new(2.1, 3.0)));
    }

    #[test]
    fn rect_empty_is_union_identity() {
        let e = Rect::empty();
        assert!(e.is_empty());
        assert_eq!(e.area(), 0.0);
        let r = Rect::new(Point::new(0.0, 0.0), Point::new(1.0, 1.0));
        assert_eq!(e.union(&r), r);
        assert!(!e.intersects(&r));
        assert!(r.contains_rect(&e));
        assert_eq!(e.min_distance_sq(Point::new(0.0, 0.0)), f64::INFINITY);
    }

    #[test]
    fn rect_intersections_are_closed() {
        let a = Rect::new(Point::new(0.0, 0.0), Point::new(1.0, 1.0));
        let b = Rect::new(Point::new(1.0, 1.0), Point::new(2.0, 2.0));
        let c = Rect::new(Point::new(1.1, 0.0), Point::new(2.0, 2.0));
        assert!(a.intersects(&b));
        assert!(!a.intersects(&c));
        assert!(a.union(&b).contains_rect(&a));
        assert_eq!(a.enlargement(&b), 3.0);
    }

    #[test]
    fn rect_min_distance() {
        let r = Rect::new(Point::new(0.0, 0.0), Point::new(1.0, 1.0));
        assert_eq!(r.min_distance_sq(Point::new(0.5, 0.5)), 0.0);
        assert_eq!(r.min_distance_sq(Point::new(3.0, 0.5)), 4.0);
        assert_eq!(r.min_distance_sq(Point::new(4.0, 5.0)), 9.0 + 16.0);
    }

    #[test]
    fn dms_and_validity() {
        let v = LonLat::dms_to_deg(-37.0, 57.0, 3.7203);
        assert!((v + 37.951_033_4).abs() < 1e-6);
        assert!(LonLat::new(3.05, 36.75).is_valid());
        assert!(!LonLat::new(181.0, 0.0).is_valid());
        assert!(!LonLat::new(0.0, f64::NAN).is_valid());
        assert_eq!(Point::from(LonLat::new(1.0, 2.0)), Point::new(1.0, 2.0));
    }
}
