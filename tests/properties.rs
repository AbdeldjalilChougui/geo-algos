//! Property tests: every fast algorithm is checked against a brute-force or
//! independent reference on random inputs.

use geo_algos::pip::{ring_contains_crossing, ring_winding_number};
use geo_algos::simplify::point_segment_distance;
use geo_algos::{
    douglas_peucker_indices, haversine, vincenty_inverse, Entry, LonLat, Location, PipAlgorithm,
    Point, Polygon, RTree, Rect, SplitStrategy,
};
use proptest::prelude::*;

fn point() -> impl Strategy<Value = Point> {
    (-100.0f64..100.0, -100.0f64..100.0).prop_map(|(x, y)| Point::new(x, y))
}

/// A random star-shaped (hence simple) polygon around the origin: one vertex
/// per angular sector, jittered within the first half of the sector, so every
/// angular gap is below pi and the origin's rays each cross exactly one edge.
fn star_polygon() -> impl Strategy<Value = Vec<Point>> {
    prop::collection::vec((0.0f64..0.5, 5.0f64..50.0), 4..40).prop_map(|v| {
        let n = v.len() as f64;
        v.iter()
            .enumerate()
            .map(|(i, &(jitter, r))| {
                let a = (i as f64 + jitter) / n * std::f64::consts::TAU;
                Point::new(a.cos() * r, a.sin() * r)
            })
            .collect()
    })
}

/// Brute-force reference for a star-shaped polygon around the origin: `p` is
/// inside iff it is closer to the origin than the boundary along its ray.
fn star_contains(ring: &[Point], p: Point) -> Option<bool> {
    let n = ring.len();
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        // Solve p = s * (a + u * (b - a)) for s > 0, u in [0, 1].
        let det = a.x * (b.y - a.y) - a.y * (b.x - a.x);
        if det.abs() < 1e-12 {
            continue;
        }
        let s = det / (p.x * (b.y - a.y) - p.y * (b.x - a.x));
        let u = (p.x * a.y - p.y * a.x) / (p.y * (b.x - a.x) - p.x * (b.y - a.y));
        if s.is_finite() && s > 0.0 && (-1e-12..=1.0 + 1e-12).contains(&u) {
            // s is the scale that maps p onto the boundary: inside if s > 1.
            if (s - 1.0).abs() < 1e-9 {
                return None; // too close to call
            }
            return Some(s > 1.0);
        }
    }
    None
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    #[test]
    fn pip_matches_star_reference(ring in star_polygon(), p in point()) {
        prop_assume!(ring.len() >= 3);
        if let Some(expected) = star_contains(&ring, p) {
            prop_assert_eq!(ring_contains_crossing(&ring, p), expected);
            prop_assert_eq!(ring_winding_number(&ring, p) != 0, expected);
        }
    }

    #[test]
    fn pip_algorithms_agree_on_simple_rings(ring in star_polygon(), p in point()) {
        prop_assume!(ring.len() >= 3);
        let poly = Polygon::without_holes(ring);
        prop_assert_eq!(
            poly.locate_with(p, PipAlgorithm::RayCasting),
            poly.locate_with(p, PipAlgorithm::WindingNumber)
        );
    }

    #[test]
    fn hole_excludes_points(p in point()) {
        let sq = |a: f64, b: f64| vec![Point::new(a, a), Point::new(b, a), Point::new(b, b), Point::new(a, b)];
        let poly = Polygon::new(sq(-50.0, 50.0), vec![sq(-10.0, 10.0)]);
        let in_outer = p.x > -50.0 && p.x < 50.0 && p.y > -50.0 && p.y < 50.0;
        let in_hole = p.x >= -10.0 && p.x <= 10.0 && p.y >= -10.0 && p.y <= 10.0;
        let loc = poly.locate(p);
        if loc != Location::Boundary {
            prop_assert_eq!(loc == Location::Inside, in_outer && !in_hole);
        }
    }

    #[test]
    fn haversine_is_a_metric(
        a in (-180.0f64..180.0, -90.0f64..90.0),
        b in (-180.0f64..180.0, -90.0f64..90.0),
        c in (-180.0f64..180.0, -90.0f64..90.0),
    ) {
        let (a, b, c) = (LonLat::new(a.0, a.1), LonLat::new(b.0, b.1), LonLat::new(c.0, c.1));
        let ab = haversine(a, b);
        prop_assert!(ab >= 0.0);
        prop_assert!((ab - haversine(b, a)).abs() < 1e-6);
        prop_assert!(ab <= haversine(a, c) + haversine(c, b) + 1e-6);
        prop_assert!(ab <= std::f64::consts::PI * 6_371_008.8 + 1e-6);
    }

    #[test]
    fn vincenty_close_to_haversine(
        a in (-180.0f64..180.0, -80.0f64..80.0),
        b in (-180.0f64..180.0, -80.0f64..80.0),
    ) {
        let (a, b) = (LonLat::new(a.0, a.1), LonLat::new(b.0, b.1));
        if let Ok(g) = vincenty_inverse(a, b) {
            let h = haversine(a, b);
            // Sphere vs ellipsoid: within 0.6% (plus a metre for tiny distances).
            prop_assert!((g.distance_m - h).abs() <= h * 0.006 + 1.0, "{} vs {}", g.distance_m, h);
        }
    }

    #[test]
    fn douglas_peucker_tolerance_holds(
        line in prop::collection::vec(point(), 0..80),
        eps in 0.0f64..20.0,
    ) {
        let idx = douglas_peucker_indices(&line, eps);
        if line.is_empty() {
            prop_assert!(idx.is_empty());
            return Ok(());
        }
        prop_assert_eq!(idx[0], 0);
        prop_assert_eq!(*idx.last().unwrap(), line.len() - 1);
        prop_assert!(idx.windows(2).all(|w| w[0] < w[1]));
        // Every dropped point is within eps of the segment that replaced it.
        for w in idx.windows(2) {
            for p in &line[w[0] + 1..w[1]] {
                prop_assert!(point_segment_distance(*p, line[w[0]], line[w[1]]) <= eps + 1e-9);
            }
        }
    }

    #[test]
    fn rtree_matches_brute_force(
        pts in prop::collection::vec(point(), 0..300),
        q1 in point(), q2 in point(), qp in point(),
        k in 1usize..12,
        bulk in any::<bool>(),
        linear in any::<bool>(),
    ) {
        let entries: Vec<Entry<usize>> = pts.iter().enumerate().map(|(i, &p)| Entry::point(p, i)).collect();
        let tree = if bulk {
            RTree::bulk_load(entries, 6)
        } else {
            let mut t = RTree::with_params(5, if linear { SplitStrategy::Linear } else { SplitStrategy::Quadratic });
            for e in entries { t.insert(e.rect, e.data); }
            t.validate(true).map_err(TestCaseError::fail)?;
            t
        };
        tree.validate(false).map_err(TestCaseError::fail)?;
        prop_assert_eq!(tree.len(), pts.len());

        // Range query.
        let q = Rect::new(q1, q2);
        let mut got: Vec<usize> = tree.search(&q).iter().map(|e| e.data).collect();
        got.sort();
        let want: Vec<usize> = (0..pts.len()).filter(|&i| q.contains_point(pts[i])).collect();
        prop_assert_eq!(got, want);

        // kNN: compare distances (ties may come in any order).
        let got: Vec<f64> = tree.nearest(qp, k).iter().map(|(d, _)| *d).collect();
        let mut all: Vec<f64> = pts.iter().map(|p| p.distance(qp)).collect();
        all.sort_by(f64::total_cmp);
        all.truncate(k);
        prop_assert_eq!(got.len(), all.len());
        for (g, w) in got.iter().zip(&all) {
            prop_assert!((g - w).abs() < 1e-9);
        }
    }
}

#[test]
fn known_distances() {
    // Algiers (Grande Poste) to Oran: ~ 351 km.
    let alger = LonLat::new(3.0588, 36.7538);
    let oran = LonLat::new(-0.6417, 35.6971);
    let h = haversine(alger, oran);
    assert!((h / 1000.0 - 351.0).abs() < 5.0, "{h}");
    // Vincenty's own test line (Flinders Peak -> Buninyong), 54972.271 m.
    let flinders = LonLat::new(
        LonLat::dms_to_deg(144.0, 25.0, 29.52440),
        LonLat::dms_to_deg(-37.0, 57.0, 3.72030),
    );
    let buninyong = LonLat::new(
        LonLat::dms_to_deg(143.0, 55.0, 35.38390),
        LonLat::dms_to_deg(-37.0, 39.0, 10.15610),
    );
    let g = vincenty_inverse(flinders, buninyong).unwrap();
    assert!((g.distance_m - 54_972.271).abs() < 0.01, "{}", g.distance_m);
}
