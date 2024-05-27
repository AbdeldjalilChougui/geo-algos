//! An R-tree (Guttman 1984) with STR bulk loading (Leutenegger et al. 1997).
//!
//! * Every node holds between `min_entries` and `max_entries` children (the
//!   root may hold fewer; STR-packed trees may have one under-full node per
//!   level, see [`RTree::bulk_load`]).
//! * All leaves are at the same depth.
//! * Each node's rectangle is the bounding box of its children.
//!
//! Nodes live in an arena (`Vec`) and refer to each other by index, which
//! keeps the borrow checker happy without `Rc<RefCell<_>>`.
//!
//! Queries:
//! * [`RTree::search`]: every entry whose rectangle intersects a query box.
//! * [`RTree::nearest`]: the `k` entries closest to a point, using best-first
//!   search (Hjaltason & Samet) with a priority queue keyed by `MINDIST`.
//!
//! Deletion is not implemented.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::types::{Point, Rect};

/// A stored item: its bounding rectangle plus a payload.
#[derive(Debug, Clone, PartialEq)]
pub struct Entry<T> {
    pub rect: Rect,
    pub data: T,
}

impl<T> Entry<T> {
    pub fn new(rect: Rect, data: T) -> Self {
        Self { rect, data }
    }

    /// Entry for a point (degenerate rectangle).
    pub fn point(p: Point, data: T) -> Self {
        Self::new(Rect::from_point(p), data)
    }
}

/// How an overflowing node is split in two.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SplitStrategy {
    /// Guttman's linear split: cheap seeds, `O(M)`.
    Linear,
    /// Guttman's quadratic split: better trees, `O(M^2)`.
    #[default]
    Quadratic,
}

#[derive(Debug, Clone)]
enum Children<T> {
    Leaf(Vec<Entry<T>>),
    Inner(Vec<usize>),
}

#[derive(Debug, Clone)]
struct Node<T> {
    rect: Rect,
    children: Children<T>,
}

impl<T> Node<T> {
    fn len(&self) -> usize {
        match &self.children {
            Children::Leaf(e) => e.len(),
            Children::Inner(c) => c.len(),
        }
    }
}

/// The R-tree.
#[derive(Debug, Clone)]
pub struct RTree<T> {
    nodes: Vec<Node<T>>,
    root: usize,
    max_entries: usize,
    min_entries: usize,
    split: SplitStrategy,
    len: usize,
}

impl<T> Default for RTree<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> RTree<T> {
    /// Empty tree with `M = 16`, `m = 6` (~40%) and quadratic split.
    pub fn new() -> Self {
        Self::with_params(16, SplitStrategy::Quadratic)
    }

    /// Empty tree with the given maximum fan-out (at least 4). The minimum
    /// fill is `40%` of it, as recommended by Guttman/Beckmann.
    pub fn with_params(max_entries: usize, split: SplitStrategy) -> Self {
        let max_entries = max_entries.max(4);
        let min_entries = ((max_entries * 2) / 5).max(2);
        Self {
            nodes: vec![Node {
                rect: Rect::empty(),
                children: Children::Leaf(Vec::new()),
            }],
            root: 0,
            max_entries,
            min_entries,
            split,
            len: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn max_entries(&self) -> usize {
        self.max_entries
    }

    /// Bounding box of everything stored (empty rect for an empty tree).
    pub fn bbox(&self) -> Rect {
        self.nodes[self.root].rect
    }

    /// Number of levels (1 for a tree that is a single leaf).
    pub fn height(&self) -> usize {
        let mut h = 1;
        let mut idx = self.root;
        while let Children::Inner(c) = &self.nodes[idx].children {
            idx = c[0];
            h += 1;
        }
        h
    }

    // ----------------------------------------------------------------- insert

    /// Inserts one entry (Guttman's ChooseLeaf / AdjustTree).
    pub fn insert(&mut self, rect: Rect, data: T) {
        // 1. Descend, remembering the path, choosing the child whose
        //    rectangle needs the least enlargement (ties: smallest area).
        let mut path = Vec::new();
        let mut idx = self.root;
        loop {
            path.push(idx);
            match &self.nodes[idx].children {
                Children::Leaf(_) => break,
                Children::Inner(children) => {
                    idx = *children
                        .iter()
                        .min_by(|&&a, &&b| {
                            let ra = &self.nodes[a].rect;
                            let rb = &self.nodes[b].rect;
                            ra.enlargement(&rect)
                                .total_cmp(&rb.enlargement(&rect))
                                .then(ra.area().total_cmp(&rb.area()))
                        })
                        .expect("inner nodes are never empty");
                }
            }
        }
        if let Children::Leaf(entries) = &mut self.nodes[idx].children {
            entries.push(Entry::new(rect, data));
        }
        self.len += 1;

        // 2. Walk back up: grow rectangles, split overflowing nodes.
        let mut split_off: Option<usize> = None;
        for depth in (0..path.len()).rev() {
            let node = path[depth];
            if let Some(sibling) = split_off.take() {
                if let Children::Inner(c) = &mut self.nodes[node].children {
                    c.push(sibling);
                }
            }
            if self.nodes[node].len() > self.max_entries {
                split_off = Some(self.split_node(node));
            } else {
                self.nodes[node].rect = self.nodes[node].rect.union(&rect);
            }
        }
        // 3. The root itself was split: grow the tree by one level.
        if let Some(sibling) = split_off {
            let old_root = self.root;
            let r = self.nodes[old_root].rect.union(&self.nodes[sibling].rect);
            self.nodes.push(Node {
                rect: r,
                children: Children::Inner(vec![old_root, sibling]),
            });
            self.root = self.nodes.len() - 1;
        }
    }

    /// Inserts a point.
    pub fn insert_point(&mut self, p: Point, data: T) {
        self.insert(Rect::from_point(p), data);
    }

    /// Splits `node` in place; returns the index of the new sibling.
    fn split_node(&mut self, node: usize) -> usize {
        let children = std::mem::replace(
            &mut self.nodes[node].children,
            Children::Inner(Vec::new()),
        );
        let (a, b) = match children {
            Children::Leaf(entries) => {
                let rects: Vec<Rect> = entries.iter().map(|e| e.rect).collect();
                let (ga, gb) = split_groups(&rects, self.min_entries, self.split);
                let (va, vb) = partition(entries, &ga, &gb);
                (Children::Leaf(va), Children::Leaf(vb))
            }
            Children::Inner(ids) => {
                let rects: Vec<Rect> = ids.iter().map(|&i| self.nodes[i].rect).collect();
                let (ga, gb) = split_groups(&rects, self.min_entries, self.split);
                let (va, vb) = partition(ids, &ga, &gb);
                (Children::Inner(va), Children::Inner(vb))
            }
        };
        self.nodes[node].children = a;
        self.nodes[node].rect = self.compute_rect(node);
        self.nodes.push(Node {
            rect: Rect::empty(),
            children: b,
        });
        let sib = self.nodes.len() - 1;
        self.nodes[sib].rect = self.compute_rect(sib);
        sib
    }

    fn compute_rect(&self, node: usize) -> Rect {
        match &self.nodes[node].children {
            Children::Leaf(e) => e.iter().fold(Rect::empty(), |r, e| r.union(&e.rect)),
            Children::Inner(c) => c
                .iter()
                .fold(Rect::empty(), |r, &i| r.union(&self.nodes[i].rect)),
        }
    }

    // -------------------------------------------------------------- bulk load

    /// Builds a packed tree with Sort-Tile-Recursive (STR).
    ///
    /// At each level the `n` items are sorted by centre `x`, cut into
    /// `ceil(sqrt(n / M))` vertical slices, each slice sorted by centre `y`
    /// and cut into nodes of `M` items. Nodes come out ~100% full, much
    /// better than repeated insertion, but the last node of each slice may be
    /// under the minimum fill.
    pub fn bulk_load(entries: Vec<Entry<T>>, max_entries: usize) -> Self {
        let mut tree = Self::with_params(max_entries, SplitStrategy::Quadratic);
        if entries.is_empty() {
            return tree;
        }
        tree.len = entries.len();
        tree.nodes.clear();
        let m = tree.max_entries;

        // Leaves.
        let groups = str_groups(entries, m, |e| e.rect.center());
        let mut level: Vec<usize> = groups
            .into_iter()
            .map(|g| {
                let rect = g.iter().fold(Rect::empty(), |r, e| r.union(&e.rect));
                tree.nodes.push(Node {
                    rect,
                    children: Children::Leaf(g),
                });
                tree.nodes.len() - 1
            })
            .collect();

        // Upper levels until a single root remains.
        while level.len() > 1 {
            let centers: Vec<(usize, Point)> = level
                .iter()
                .map(|&i| (i, tree.nodes[i].rect.center()))
                .collect();
            let groups = str_groups(centers, m, |&(_, c)| c);
            level = groups
                .into_iter()
                .map(|g| {
                    let ids: Vec<usize> = g.into_iter().map(|(i, _)| i).collect();
                    let rect = ids
                        .iter()
                        .fold(Rect::empty(), |r, &i| r.union(&tree.nodes[i].rect));
                    tree.nodes.push(Node {
                        rect,
                        children: Children::Inner(ids),
                    });
                    tree.nodes.len() - 1
                })
                .collect();
        }
        tree.root = level[0];
        tree
    }

    // ---------------------------------------------------------------- queries

    /// All entries whose rectangle intersects `query` (closed boundaries).
    pub fn search(&self, query: &Rect) -> Vec<&Entry<T>> {
        let mut out = Vec::new();
        if self.is_empty() {
            return out;
        }
        let mut stack = vec![self.root];
        while let Some(idx) = stack.pop() {
            let node = &self.nodes[idx];
            if !node.rect.intersects(query) {
                continue;
            }
            match &node.children {
                Children::Leaf(entries) => {
                    out.extend(entries.iter().filter(|e| e.rect.intersects(query)))
                }
                Children::Inner(c) => stack.extend(c.iter().copied()),
            }
        }
        out
    }

    /// The `k` entries nearest to `p`, closest first, with their Euclidean
    /// distance (distance to the entry's rectangle, 0 if `p` is inside it).
    ///
    /// Best-first search: a min-heap mixes nodes and entries, keyed by
    /// `MINDIST`. Because a node's `MINDIST` is a lower bound for everything
    /// below it, an entry popped from the heap is guaranteed to be the next
    /// nearest one.
    pub fn nearest(&self, p: Point, k: usize) -> Vec<(f64, &Entry<T>)> {
        let mut out = Vec::with_capacity(k);
        if k == 0 || self.is_empty() {
            return out;
        }
        let mut heap = BinaryHeap::new();
        heap.push(HeapItem {
            dist_sq: self.nodes[self.root].rect.min_distance_sq(p),
            item: Item::Node(self.root),
        });
        while let Some(HeapItem { dist_sq, item }) = heap.pop() {
            match item {
                Item::Entry(node, i) => {
                    if let Children::Leaf(entries) = &self.nodes[node].children {
                        out.push((dist_sq.sqrt(), &entries[i]));
                    }
                    if out.len() == k {
                        break;
                    }
                }
                Item::Node(idx) => match &self.nodes[idx].children {
                    Children::Leaf(entries) => {
                        for (i, e) in entries.iter().enumerate() {
                            heap.push(HeapItem {
                                dist_sq: e.rect.min_distance_sq(p),
                                item: Item::Entry(idx, i),
                            });
                        }
                    }
                    Children::Inner(c) => {
                        for &child in c {
                            heap.push(HeapItem {
                                dist_sq: self.nodes[child].rect.min_distance_sq(p),
                                item: Item::Node(child),
                            });
                        }
                    }
                },
            }
        }
        out
    }

    /// Iterates over every stored entry (no particular order).
    pub fn iter(&self) -> impl Iterator<Item = &Entry<T>> {
        self.nodes.iter().flat_map(|n| match &n.children {
            Children::Leaf(e) => e.iter(),
            Children::Inner(_) => [].iter(),
        })
    }

    /// Checks the structural invariants; returns a description of the first
    /// violation. `strict_min_fill` also enforces the minimum fill on
    /// non-root nodes (true for insert-built trees, not for STR).
    pub fn validate(&self, strict_min_fill: bool) -> Result<(), String> {
        let mut leaf_depth = None;
        let mut count = 0;
        let mut stack = vec![(self.root, 0usize)];
        while let Some((idx, depth)) = stack.pop() {
            let node = &self.nodes[idx];
            let n = node.len();
            if n > self.max_entries {
                return Err(format!("node {idx} has {n} > M children"));
            }
            if idx != self.root && n == 0 {
                return Err(format!("node {idx} is empty"));
            }
            if strict_min_fill && idx != self.root && n < self.min_entries {
                return Err(format!("node {idx} has {n} < m children"));
            }
            if node.rect != self.compute_rect(idx) {
                return Err(format!("node {idx} rect is not the bbox of its children"));
            }
            match &node.children {
                Children::Leaf(e) => {
                    count += e.len();
                    match leaf_depth {
                        None => leaf_depth = Some(depth),
                        Some(d) if d != depth => {
                            return Err(format!("leaves at depths {d} and {depth}"))
                        }
                        _ => {}
                    }
                }
                Children::Inner(c) => stack.extend(c.iter().map(|&i| (i, depth + 1))),
            }
        }
        if count != self.len {
            return Err(format!("{count} entries reachable, len is {}", self.len));
        }
        Ok(())
    }
}

// --------------------------------------------------------------------- helpers

#[derive(Debug, Clone, Copy, PartialEq)]
enum Item {
    Node(usize),
    Entry(usize, usize),
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct HeapItem {
    dist_sq: f64,
    item: Item,
}

impl Eq for HeapItem {}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reversed: BinaryHeap is a max-heap, we want the smallest distance.
        // On ties, entries come out before nodes so results are produced early.
        other.dist_sq.total_cmp(&self.dist_sq).then_with(|| {
            let rank = |i: &Item| matches!(i, Item::Entry(..)) as u8;
            rank(&self.item).cmp(&rank(&other.item))
        })
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Moves `items` into two vectors following the index groups.
fn partition<X>(items: Vec<X>, ga: &[usize], gb: &[usize]) -> (Vec<X>, Vec<X>) {
    let mut slots: Vec<Option<X>> = items.into_iter().map(Some).collect();
    let a = ga.iter().map(|&i| slots[i].take().unwrap()).collect();
    let b = gb.iter().map(|&i| slots[i].take().unwrap()).collect();
    (a, b)
}

/// Splits `rects` (indices) into two groups, each with at least `min` items.
fn split_groups(rects: &[Rect], min: usize, strategy: SplitStrategy) -> (Vec<usize>, Vec<usize>) {
    let (s1, s2) = match strategy {
        SplitStrategy::Quadratic => quadratic_seeds(rects),
        SplitStrategy::Linear => linear_seeds(rects),
    };
    let mut ga = vec![s1];
    let mut gb = vec![s2];
    let mut ra = rects[s1];
    let mut rb = rects[s2];
    let mut rest: Vec<usize> = (0..rects.len()).filter(|&i| i != s1 && i != s2).collect();

    while !rest.is_empty() {
        // If one group needs every remaining item to reach the minimum, give them all.
        if ga.len() + rest.len() == min {
            ga.append(&mut rest);
            break;
        }
        if gb.len() + rest.len() == min {
            gb.append(&mut rest);
            break;
        }
        // PickNext: quadratic picks the item with the strongest preference,
        // linear just takes the next one.
        let pos = match strategy {
            SplitStrategy::Quadratic => (0..rest.len())
                .max_by(|&x, &y| {
                    let d = |i: usize| {
                        let r = &rects[rest[i]];
                        (ra.enlargement(r) - rb.enlargement(r)).abs()
                    };
                    d(x).total_cmp(&d(y))
                })
                .unwrap(),
            SplitStrategy::Linear => rest.len() - 1,
        };
        let i = rest.swap_remove(pos);
        let r = &rects[i];
        let (ea, eb) = (ra.enlargement(r), rb.enlargement(r));
        let to_a = match ea.total_cmp(&eb) {
            Ordering::Less => true,
            Ordering::Greater => false,
            Ordering::Equal => match ra.area().total_cmp(&rb.area()) {
                Ordering::Less => true,
                Ordering::Greater => false,
                Ordering::Equal => ga.len() <= gb.len(),
            },
        };
        if to_a {
            ga.push(i);
            ra = ra.union(r);
        } else {
            gb.push(i);
            rb = rb.union(r);
        }
    }
    (ga, gb)
}

/// Quadratic PickSeeds: the pair that would waste the most area together.
fn quadratic_seeds(rects: &[Rect]) -> (usize, usize) {
    let mut best = (0, 1);
    let mut worst = f64::NEG_INFINITY;
    for i in 0..rects.len() {
        for j in i + 1..rects.len() {
            let d = rects[i].union(&rects[j]).area() - rects[i].area() - rects[j].area();
            if d > worst {
                worst = d;
                best = (i, j);
            }
        }
    }
    best
}

/// Linear PickSeeds: along each axis, the rectangle with the highest low side
/// and the one with the lowest high side; keep the axis with the largest
/// separation normalised by the total extent.
fn linear_seeds(rects: &[Rect]) -> (usize, usize) {
    let mut best = (0, 1);
    let mut best_sep = f64::NEG_INFINITY;
    for axis in 0..2 {
        let lo = |r: &Rect| if axis == 0 { r.min.x } else { r.min.y };
        let hi = |r: &Rect| if axis == 0 { r.max.x } else { r.max.y };
        let highest_low = (0..rects.len())
            .max_by(|&a, &b| lo(&rects[a]).total_cmp(&lo(&rects[b])))
            .unwrap();
        let lowest_high = (0..rects.len())
            .filter(|&i| i != highest_low)
            .min_by(|&a, &b| hi(&rects[a]).total_cmp(&hi(&rects[b])))
            .unwrap();
        let min_lo = rects.iter().map(lo).fold(f64::INFINITY, f64::min);
        let max_hi = rects.iter().map(hi).fold(f64::NEG_INFINITY, f64::max);
        let extent = (max_hi - min_lo).max(f64::MIN_POSITIVE);
        let sep = (lo(&rects[highest_low]) - hi(&rects[lowest_high])) / extent;
        if sep > best_sep {
            best_sep = sep;
            best = (lowest_high, highest_low);
        }
    }
    best
}

/// Sort-Tile-Recursive grouping of `items` into chunks of at most `m`.
fn str_groups<X>(mut items: Vec<X>, m: usize, center: impl Fn(&X) -> Point) -> Vec<Vec<X>> {
    let n = items.len();
    let pages = n.div_ceil(m);
    let slices = (pages as f64).sqrt().ceil() as usize;
    let slice_size = slices * m;
    items.sort_by(|a, b| center(a).x.total_cmp(&center(b).x));
    let mut groups = Vec::with_capacity(pages);
    let mut iter = items.into_iter().peekable();
    while iter.peek().is_some() {
        let mut slice: Vec<X> = iter.by_ref().take(slice_size).collect();
        slice.sort_by(|a, b| center(a).y.total_cmp(&center(b).y));
        let mut s = slice.into_iter().peekable();
        while s.peek().is_some() {
            groups.push(s.by_ref().take(m).collect());
        }
    }
    groups
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(n: usize) -> Vec<Entry<usize>> {
        (0..n * n)
            .map(|i| Entry::point(Point::new((i % n) as f64, (i / n) as f64), i))
            .collect()
    }

    #[test]
    fn empty_tree() {
        let t: RTree<()> = RTree::new();
        assert!(t.is_empty());
        assert!(t.search(&Rect::new(Point::new(0.0, 0.0), Point::new(1.0, 1.0))).is_empty());
        assert!(t.nearest(Point::new(0.0, 0.0), 3).is_empty());
        assert_eq!(t.height(), 1);
        t.validate(true).unwrap();
    }

    #[test]
    fn insert_grows_and_stays_valid() {
        for split in [SplitStrategy::Linear, SplitStrategy::Quadratic] {
            let mut t = RTree::with_params(4, split);
            for e in grid(20) {
                t.insert(e.rect, e.data);
                t.validate(true).unwrap();
            }
            assert_eq!(t.len(), 400);
            assert!(t.height() >= 4);
            assert_eq!(t.iter().count(), 400);
        }
    }

    #[test]
    fn range_query() {
        let t = RTree::bulk_load(grid(10), 8);
        t.validate(false).unwrap();
        let q = Rect::new(Point::new(2.0, 3.0), Point::new(4.0, 4.0));
        let mut found: Vec<usize> = t.search(&q).into_iter().map(|e| e.data).collect();
        found.sort();
        assert_eq!(found, vec![32, 33, 34, 42, 43, 44]);
    }

    #[test]
    fn knn_query() {
        let t = RTree::bulk_load(grid(10), 4);
        let res = t.nearest(Point::new(5.1, 5.0), 3);
        assert_eq!(res[0].1.data, 55);
        assert!((res[0].0 - 0.1).abs() < 1e-12);
        assert_eq!(res.len(), 3);
        assert!(res.windows(2).all(|w| w[0].0 <= w[1].0));
        assert_eq!(t.nearest(Point::new(0.0, 0.0), 1000).len(), 100);
    }

    #[test]
    fn rect_entries() {
        let mut t = RTree::new();
        t.insert(Rect::new(Point::new(0.0, 0.0), Point::new(10.0, 10.0)), "big");
        t.insert(Rect::new(Point::new(20.0, 20.0), Point::new(21.0, 21.0)), "small");
        let hit = t.search(&Rect::from_point(Point::new(5.0, 5.0)));
        assert_eq!(hit.len(), 1);
        assert_eq!(hit[0].data, "big");
        let nn = t.nearest(Point::new(5.0, 5.0), 2);
        assert_eq!(nn[0].0, 0.0);
        assert_eq!(nn[1].1.data, "small");
    }
}
