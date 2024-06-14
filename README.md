# geo-algos

Classic geospatial algorithms written from scratch in Rust, with no runtime
dependencies, for learning. Each one is checked against a brute-force or
independent reference with property tests.

| Module | What it does |
|---|---|
| `pip` | Point in polygon with holes: ray casting (even-odd) and winding number (non-zero), explicit `Boundary` result |
| `distance` | Haversine (sphere) and Vincenty inverse (WGS84 ellipsoid) with bearings and a convergence error |
| `simplify` | Douglas-Peucker polyline simplification, iterative (no recursion depth limit) |
| `rtree` | R-tree: Guttman insert with linear or quadratic split, STR bulk load, range query, best-first kNN |

## Install

```sh
cargo add geo-algos
```

## Quick start

```sh
cargo test                      # unit tests + proptest checks
cargo run --example demo        # small tour of every module
cargo bench                     # criterion benchmarks
```

```rust
use geo_algos::{haversine, LonLat, Point, Polygon, RTree, Entry, Rect};

let d = haversine(LonLat::new(3.0588, 36.7538), LonLat::new(-0.6417, 35.6971)); // ~352 km

let sq = |a: f64, b: f64| vec![Point::new(a, a), Point::new(b, a), Point::new(b, b), Point::new(a, b)];
let donut = Polygon::new(sq(0.0, 10.0), vec![sq(4.0, 6.0)]);
assert!(donut.contains(Point::new(1.0, 1.0)));
assert!(!donut.contains(Point::new(5.0, 5.0))); // in the hole

let tree = RTree::bulk_load(vec![Entry::point(Point::new(1.0, 2.0), "a")], 16);
let near = tree.nearest(Point::new(0.0, 0.0), 1);
```

## How it works

- **Point in polygon.** Ray casting counts edge crossings of a horizontal ray
  with a half-open rule, so a ray through a vertex is counted once. The
  winding number (Sunday's version) uses only orientation tests. Both agree on
  simple rings and differ on self-intersecting ones (the centre of a pentagram
  is outside under even-odd and inside under non-zero; there is a test for
  it). `Polygon::locate` checks the rings first, so points exactly on an edge
  come back as `Boundary` instead of an arbitrary answer.
- **Distances.** Haversine on a sphere of mean radius 6,371,008.8 m is off by
  up to about 0.5% compared with the ellipsoid. Vincenty iterates on the
  longitude difference on WGS84. It gives sub-millimetre accuracy and is
  checked against Vincenty's own Flinders Peak to Buninyong line
  (54,972.271 m). Nearly antipodal points can fail to converge, which is
  returned as `VincentyError::NoConvergence`.
- **Douglas-Peucker.** Keeps the endpoints, finds the farthest point from the
  chord and recurses while it is beyond `epsilon`. An explicit stack replaces
  recursion. The property test checks the guarantee: every dropped point is
  within `epsilon` of the segment that replaced it.
- **R-tree.** Nodes are stored in an arena (`Vec`) and referenced by index.
  `insert` follows Guttman: ChooseLeaf by least enlargement, then split on
  overflow (quadratic PickSeeds/PickNext or linear seeds) and adjust up to the
  root. `bulk_load` is Sort-Tile-Recursive, which packs nodes full. `nearest`
  is best-first search with a min-heap keyed by MINDIST, so entries come out
  in exact distance order. `validate()` checks the invariants: equal leaf
  depth, fan-out bounds and tight bounding boxes.

## Tests

- 39 unit tests next to the code: edge cases, boundaries, holes, a pentagram,
  degenerate inputs and known distances.
- 8 proptest properties in `tests/properties.rs`:
  - PIP against an analytic star-polygon reference, plus ray casting vs winding number
  - holes against a closed-form predicate
  - Haversine metric axioms
  - Vincenty within 0.6% of Haversine
  - the Douglas-Peucker tolerance guarantee
  - R-tree range and kNN results against brute force, for insert-built trees (both splits) and STR trees

## Benchmarks

`cargo bench -- --quick` on a laptop that was busy with other builds, so treat
these as orders of magnitude only:

| benchmark | time |
|---|---|
| PIP, 1000-vertex ring | ~58 µs |
| Haversine | ~0.23 µs |
| Vincenty | ~2 µs |
| Douglas-Peucker, 10k points | ~12 ms |
| R-tree STR bulk load, 50k points | ~52 ms |
| R-tree one-by-one insert, 50k points | ~0.3–0.5 s |
| R-tree range query (50×50 box over 50k points) | ~30 µs |
| R-tree 10-NN | ~0.1 ms |

## Status and limitations

- This is a learning project. For real work use the `geo` and `rstar` crates.
- Everything is plain `f64` with no robust predicates. Boundary detection is
  exact only when the arithmetic is: vertices and axis-aligned edges work, but
  a point interpolated along a slanted edge may land a hair inside or outside.
- Planar algorithms treat lon/lat degrees as flat x/y. Nothing handles the
  antimeridian, and Douglas-Peucker `epsilon` is in input units.
- The R-tree has no deletion and no R*-tree reinsertion. kNN distance is
  measured to each entry's bounding rectangle.
- Hole validity (holes inside the exterior, not overlapping) is assumed and
  not checked.

## License

MIT, see [LICENSE](LICENSE). No external data is used.
