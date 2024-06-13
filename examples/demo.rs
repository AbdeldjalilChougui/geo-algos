//! `cargo run --example demo`
use geo_algos::{
    douglas_peucker, haversine, vincenty_inverse, Entry, LonLat, Point, Polygon, RTree, Rect,
};

fn main() {
    let alger = LonLat::new(3.0588, 36.7538);
    let oran = LonLat::new(-0.6417, 35.6971);
    println!("Alger -> Oran haversine: {:.1} km", haversine(alger, oran) / 1000.0);
    let g = vincenty_inverse(alger, oran).expect("converges");
    println!(
        "Alger -> Oran vincenty:  {:.1} km, initial bearing {:.1} deg ({} iterations)",
        g.distance_m / 1000.0,
        g.initial_bearing_deg,
        g.iterations
    );

    let sq = |a: f64, b: f64| vec![Point::new(a, a), Point::new(b, a), Point::new(b, b), Point::new(a, b)];
    let donut = Polygon::new(sq(0.0, 10.0), vec![sq(4.0, 6.0)]);
    for p in [Point::new(1.0, 1.0), Point::new(5.0, 5.0), Point::new(4.0, 5.0)] {
        println!("{p:?} -> {:?}", donut.locate(p));
    }

    let wiggly: Vec<Point> = (0..100).map(|i| Point::new(i as f64, (i as f64 / 5.0).sin())).collect();
    println!("Douglas-Peucker: {} -> {} points (eps 0.1)", wiggly.len(), douglas_peucker(&wiggly, 0.1).len());

    let tree = RTree::bulk_load(
        (0..10_000).map(|i| Entry::point(Point::new((i % 100) as f64, (i / 100) as f64), i)).collect(),
        16,
    );
    let hits = tree.search(&Rect::new(Point::new(10.0, 10.0), Point::new(12.0, 11.0)));
    println!("R-tree: height {}, {} hits in range", tree.height(), hits.len());
    for (d, e) in tree.nearest(Point::new(50.2, 50.1), 3) {
        println!("  nearest #{} at {:.3}", e.data, d);
    }
}
