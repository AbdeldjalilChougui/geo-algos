//! Douglas-Peucker (Ramer-Douglas-Peucker) polyline simplification.
//!
//! Given a polyline and a tolerance `epsilon`, keep the two endpoints, find the
//! point farthest from the segment joining them and, if it is farther than
//! `epsilon`, keep it and recurse on both halves. Otherwise drop everything in
//! between.
//!
//! Guarantee: every removed point lies within `epsilon` of the segment of the
//! simplified line that replaced it. The output is a subsequence of the input
//! and always keeps the first and last points.
//!
//! The recursion is done with an explicit stack, so long lines cannot overflow
//! the call stack. Worst case is `O(n^2)`, typical case `O(n log n)`.
//!
//! Distances are planar (Euclidean on `x`/`y`). For lon/lat input, `epsilon`
//! is therefore in degrees.

use crate::types::Point;

/// Distance from `p` to the closed segment `[a, b]`.
pub fn point_segment_distance(p: Point, a: Point, b: Point) -> f64 {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let len_sq = dx * dx + dy * dy;
    if len_sq == 0.0 {
        return p.distance(a);
    }
    let t = (((p.x - a.x) * dx + (p.y - a.y) * dy) / len_sq).clamp(0.0, 1.0);
    p.distance(Point::new(a.x + t * dx, a.y + t * dy))
}

/// Indices of the points kept by Douglas-Peucker, in increasing order.
///
/// Lines with fewer than 3 points are returned unchanged. A negative or NaN
/// `epsilon` is treated as `0.0` (only exactly-collinear points are removed).
pub fn douglas_peucker_indices(points: &[Point], epsilon: f64) -> Vec<usize> {
    let n = points.len();
    if n < 3 {
        return (0..n).collect();
    }
    let eps = if epsilon.is_nan() { 0.0 } else { epsilon.max(0.0) };
    let mut keep = vec![false; n];
    keep[0] = true;
    keep[n - 1] = true;

    let mut stack = vec![(0usize, n - 1)];
    while let Some((start, end)) = stack.pop() {
        if end <= start + 1 {
            continue;
        }
        let (a, b) = (points[start], points[end]);
        let mut max_d = -1.0;
        let mut max_i = start;
        for (i, &p) in points.iter().enumerate().take(end).skip(start + 1) {
            let d = point_segment_distance(p, a, b);
            if d > max_d {
                max_d = d;
                max_i = i;
            }
        }
        if max_d > eps {
            keep[max_i] = true;
            stack.push((start, max_i));
            stack.push((max_i, end));
        }
    }
    keep.iter()
        .enumerate()
        .filter_map(|(i, &k)| k.then_some(i))
        .collect()
}

/// Simplified copy of the polyline (see [`douglas_peucker_indices`]).
pub fn douglas_peucker(points: &[Point], epsilon: f64) -> Vec<Point> {
    douglas_peucker_indices(points, epsilon)
        .into_iter()
        .map(|i| points[i])
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pts(v: &[(f64, f64)]) -> Vec<Point> {
        v.iter().map(|&(x, y)| Point::new(x, y)).collect()
    }

    #[test]
    fn segment_distance() {
        let a = Point::new(0.0, 0.0);
        let b = Point::new(10.0, 0.0);
        assert_eq!(point_segment_distance(Point::new(5.0, 3.0), a, b), 3.0);
        // Beyond the endpoint: distance to the endpoint, not to the infinite line.
        assert_eq!(point_segment_distance(Point::new(13.0, 4.0), a, b), 5.0);
        assert_eq!(point_segment_distance(Point::new(3.0, 4.0), a, a), 5.0);
    }

    #[test]
    fn short_lines_unchanged() {
        assert!(douglas_peucker(&[], 1.0).is_empty());
        let two = pts(&[(0.0, 0.0), (1.0, 1.0)]);
        assert_eq!(douglas_peucker(&two, 1.0), two);
    }

    #[test]
    fn collinear_points_removed() {
        let line = pts(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]);
        assert_eq!(douglas_peucker_indices(&line, 0.0), vec![0, 3]);
    }

    #[test]
    fn classic_example() {
        let line = pts(&[
            (0.0, 0.0),
            (1.0, 0.1),
            (2.0, -0.1),
            (3.0, 5.0),
            (4.0, 6.0),
            (5.0, 7.0),
            (6.0, 8.1),
            (7.0, 9.0),
            (8.0, 9.0),
            (9.0, 9.0),
        ]);
        let idx = douglas_peucker_indices(&line, 1.0);
        assert_eq!(idx.first(), Some(&0));
        assert_eq!(idx.last(), Some(&9));
        assert!(idx.contains(&2) || idx.contains(&3));
        assert!(idx.len() < line.len());
        // Huge epsilon keeps only endpoints, zero keeps every non-collinear point.
        assert_eq!(douglas_peucker_indices(&line, 100.0), vec![0, 9]);
        assert!(!douglas_peucker_indices(&line, 0.0).contains(&8)); // (8,9) is collinear
    }

    #[test]
    fn negative_or_nan_epsilon() {
        let line = pts(&[(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)]);
        assert_eq!(douglas_peucker(&line, -1.0).len(), 3);
        assert_eq!(douglas_peucker(&line, f64::NAN).len(), 3);
    }

    #[test]
    fn long_line_does_not_overflow() {
        let line: Vec<Point> = (0..5_000)
            .map(|i| Point::new(i as f64, if i % 2 == 0 { 0.0 } else { 1.0 }))
            .collect();
        assert_eq!(douglas_peucker(&line, 0.5).len(), line.len());
        assert_eq!(douglas_peucker(&line, 2.0).len(), 2);
    }
}
