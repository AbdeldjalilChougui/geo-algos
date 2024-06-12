//! Criterion benchmarks: `cargo bench`.
use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion};
use geo_algos::{douglas_peucker, haversine, vincenty_inverse, Entry, LonLat, Point, Polygon, RTree, Rect};

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    (*seed >> 11) as f64 / (1u64 << 53) as f64
}

fn bench(c: &mut Criterion) {
    let mut s = 42;
    let circle: Vec<Point> = (0..1000)
        .map(|i| {
            let a = i as f64 / 1000.0 * std::f64::consts::TAU;
            Point::new(a.cos() * 100.0, a.sin() * 100.0)
        })
        .collect();
    let poly = Polygon::without_holes(circle);
    c.bench_function("pip winding 1000-vertex ring", |b| {
        b.iter(|| poly.locate(black_box(Point::new(10.0, 20.0))))
    });

    let (a, z) = (LonLat::new(3.0588, 36.7538), LonLat::new(-0.6417, 35.6971));
    c.bench_function("haversine", |b| b.iter(|| haversine(black_box(a), black_box(z))));
    c.bench_function("vincenty", |b| b.iter(|| vincenty_inverse(black_box(a), black_box(z))));

    let line: Vec<Point> = (0..10_000)
        .map(|i| Point::new(i as f64, (i as f64 / 50.0).sin() * 10.0 + lcg(&mut s)))
        .collect();
    c.bench_function("douglas-peucker 10k points", |b| b.iter(|| douglas_peucker(black_box(&line), 0.5)));

    let pts: Vec<Point> = (0..50_000).map(|_| Point::new(lcg(&mut s) * 1000.0, lcg(&mut s) * 1000.0)).collect();
    let entries: Vec<Entry<usize>> = pts.iter().enumerate().map(|(i, &p)| Entry::point(p, i)).collect();
    c.bench_function("rtree STR bulk load 50k", |b| b.iter(|| RTree::bulk_load(entries.clone(), 16)));
    c.bench_function("rtree insert 50k", |b| {
        b.iter(|| {
            let mut t = RTree::new();
            for (i, &p) in pts.iter().enumerate() {
                t.insert_point(p, i);
            }
            t
        })
    });
    let t = RTree::bulk_load(entries, 16);
    let q = Rect::new(Point::new(100.0, 100.0), Point::new(150.0, 150.0));
    c.bench_function("rtree range query", |b| b.iter(|| t.search(black_box(&q)).len()));
    c.bench_function("rtree 10-NN", |b| b.iter(|| t.nearest(black_box(Point::new(500.0, 500.0)), 10)));
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(20);
    targets = bench
}
criterion_main!(benches);
