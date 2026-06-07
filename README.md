# graph-homology

> **The topology hidden in graphs. Clique complexes, graph Laplacians, Euler characteristic.**

[![crates.io](https://img.shields.io/crates/v/graph-homology.svg)](https://crates.io/crates/graph-homology)
[![docs.rs](https://docs.rs/graph-homology/badge.svg)](https://docs.rs/graph-homology)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A Rust library computing homological invariants of graphs and their clique complexes. Constructs simplicial complexes from graph cliques, computes Betti numbers, Euler characteristic, and graph Laplacians (combinatorial and normalized). Turns any graph into a topological space and extracts its shape.

---

## Table of Contents

- [What is Graph Homology?](#what-is-graph-homology)
- [Why Does This Matter?](#why-does-this-matter)
- [Architecture](#architecture)
- [Quick Start](#quick-start)
- [API Reference](#api-reference)
- [Mathematical Background](#mathematical-background)
- [Installation](#installation)
- [Related Crates](#related-crates)
- [License](#license)

---

## What is Graph Homology?

Every graph has a hidden topological structure. The **clique complex** Cl(G) of a graph G is the simplicial complex whose simplices are the complete subgraphs (cliques) of G:

```
Graph G:           Clique Complex Cl(G):
●──●──●            Vertices:  0-simplices {a}, {b}, {c}, {d}
│  │  │            Edges:     1-simplices {a,b}, {b,c}, {c,d}, {a,b,c}
●──●──●            Triangles: 2-simplex  {a,b,c}
```

From this complex, we extract **homological invariants**:

- **H₀** (Betti number β₀) = number of connected components
- **H₁** (Betti number β₁) = number of independent cycles (holes not filled by triangles)
- **Euler characteristic** χ = V − E + F = Σ(−1)ᵏ fₖ
- **Graph Laplacian** eigenvalues → spectral structure

```
Graph:    ●──●──●           β₀ = 1, β₁ = 0, χ = 3−2 = 1
          │  │
          ●──●

Graph:    ●──●              β₀ = 1, β₁ = 1, χ = 4−5+1 = 0
          │  │              (one cycle: the square)
          ●──●

Graph:    ●──●──●           β₀ = 2, β₁ = 0, χ = 6−5 = 1
          │     │           (two disconnected components)
          ●──●──●
```

## Why Does This Matter?

**For network analysis**: Betti numbers reveal the "shape" of a network — how many separate clusters (β₀) and how many loop structures (β₁) exist. This goes beyond degree distributions and clustering coefficients.

**For graph neural networks**: The combinatorial Laplacian is the foundation of spectral GNNs. Normalized Laplacian eigenvalues determine mixing time, clustering quality, and information diffusion speed.

**For multi-agent systems**: The topology of agent communication graphs determines consensus reachability, fault tolerance, and information propagation. A disconnected graph (β₀ > 1) cannot reach consensus.

**For topological data analysis**: Graph homology is the starting point for persistent homology on network data — tracking how topological features appear and disappear as edges are added.

## Architecture

```
graph-homology
│
├── SimplicialGraph            ← Undirected graph representation
│   ├── new() / with_vertices(n)   Empty / n-vertex graph
│   ├── add_vertex(v)              Add vertex
│   ├── add_edge(u, v)             Add undirected edge
│   ├── neighbors(v)               Adjacent vertices
│   ├── complete(n)                Complete graph K_n
│   ├── cycle(n)                   Cycle graph C_n
│   ├── path(n)                    Path graph P_n
│   ├── connected_components()     Find all components
│   └── degree(v)                  Vertex degree
│
├── CliqueComplex              ← Simplicial complex from cliques
│   ├── from_graph(graph)          Build Cl(G)
│   ├── simplices(k)               All k-simplices
│   ├── simplex_count(k)           Number of k-simplices
│   └── max_dimension()            Highest simplex dimension
│
├── GraphHomology              ← Homology computation
│   ├── h0(graph)                  Connected components
│   ├── h1(graph)                  Independent cycles
│   └── betti_numbers(complex, k)  β₀, β₁, ..., βₖ
│
├── EulerCharacteristic        ← Euler characteristic
│   ├── compute(graph)             χ = V − E + F
│   ├── compute_complex(complex)   χ = Σ(−1)ᵏ fₖ
│   └── count_triangles(graph)     Number of 3-cliques
│
└── GraphLaplacian             ← Spectral graph theory
    ├── combinatorial(graph)       L = D − A
    ├── normalized(graph)          L_norm = D^{-1/2} L D^{-1/2}
    ├── eigenvalues_approx(L, n)   Power iteration eigenvalues
    ├── zero_eigenvalue_count(graph) = β₀
    └── algebraic_connectivity(graph)  Second-smallest eigenvalue
```

## Quick Start

```rust
use graph_homology::{
    SimplicialGraph, CliqueComplex, GraphHomology,
    EulerCharacteristic, GraphLaplacian,
};

// Build a square graph: 4 nodes, 4 edges, 0 triangles
let mut g = SimplicialGraph::new();
g.add_edge(0, 1);
g.add_edge(1, 2);
g.add_edge(2, 3);
g.add_edge(3, 0);

// Basic graph properties
println!("Vertices: {}", g.vertex_count());
println!("Edges: {}", g.edge_count());
println!("Components: {:?}", g.connected_components());

// Compute homology
let h0 = GraphHomology::h0(&g);  // 1 (one component)
let h1 = GraphHomology::h1(&g);  // 1 (one cycle)
println!("β₀ = {}, β₁ = {}", h0, h1);

// Build clique complex
let complex = CliqueComplex::from_graph(g.clone());
println!("Max dimension: {}", complex.max_dimension());
println!("Triangles: {}", complex.simplex_count(2));

// Euler characteristic: χ = V − E + F
let chi = EulerCharacteristic::compute(&g);
println!("Euler characteristic: {}", chi); // 4 - 4 + 0 = 0

// Graph Laplacian: L = D - A
let L = GraphLaplacian::combinatorial(&g);
let L_norm = GraphLaplacian::normalized(&g);
println!("Zero eigenvalues: {}", GraphLaplacian::zero_eigenvalue_count(&g));
println!("Algebraic connectivity: {:.4}", GraphLaplacian::algebraic_connectivity(&g));

// Pre-built graphs
let K5 = SimplicialGraph::complete(5);
let C6 = SimplicialGraph::cycle(6);
let P4 = SimplicialGraph::path(4);
println!("K₅ β₁ = {}", GraphHomology::h1(&K5)); // 0 (complete = no holes)
println!("C₆ β₁ = {}", GraphHomology::h1(&C6)); // 1 (one cycle)
```

## API Reference

### SimplicialGraph

| Method | Returns | Description |
|--------|---------|-------------|
| `new()` | `Self` | Empty graph |
| `with_vertices(n)` | `Self` | n isolated vertices |
| `add_vertex(v)` | `()` | Add vertex |
| `add_edge(u, v)` | `()` | Add undirected edge |
| `neighbors(v)` | `Option<&HashSet<usize>>` | Adjacent vertices |
| `vertices()` | `Vec<usize>` | All vertices |
| `edges()` | `Vec<(usize, usize)>` | All edges |
| `vertex_count()` | `usize` | Number of vertices |
| `edge_count()` | `usize` | Number of edges |
| `has_edge(u, v)` | `bool` | Edge existence |
| `degree(v)` | `usize` | Vertex degree |
| `complete(n)` | `Self` | Complete graph K_n |
| `cycle(n)` | `Self` | Cycle graph C_n |
| `path(n)` | `Self` | Path graph P_n |
| `connected_components()` | `Vec<Vec<usize>>` | All components |

### CliqueComplex

| Method | Returns | Description |
|--------|---------|-------------|
| `from_graph(graph)` | `Self` | Build from graph |
| `simplices(k)` | `Vec<BTreeSet<usize>>` | All k-simplices |
| `simplex_count(k)` | `usize` | Number of k-simplices |
| `max_dimension()` | `usize` | Highest dimension |

### GraphHomology

| Method | Returns | Description |
|--------|---------|-------------|
| `h0(graph)` | `usize` | Number of connected components (β₀) |
| `h1(graph)` | `usize` | Number of independent cycles (β₁) |
| `betti_numbers(complex, max_dim)` | `Vec<usize>` | β₀, β₁, ..., βₖ |

### EulerCharacteristic

| Method | Returns | Description |
|--------|---------|-------------|
| `compute(graph)` | `i32` | χ = V − E + F |
| `compute_complex(complex)` | `i32` | χ = Σ(−1)ᵏ fₖ |
| `count_triangles(graph)` | `usize` | Number of 3-cliques |

### GraphLaplacian

| Method | Returns | Description |
|--------|---------|-------------|
| `combinatorial(graph)` | `Vec<Vec<f64>>` | L = D − A |
| `normalized(graph)` | `Vec<Vec<f64>>` | D^{-1/2} L D^{-1/2} |
| `eigenvalues_approx(L, iters)` | `Vec<f64>` | Power iteration |
| `zero_eigenvalue_count(graph)` | `usize` | = β₀ (components) |
| `algebraic_connectivity(graph)` | `f64` | 2nd-smallest eigenvalue |

## Mathematical Background

### Clique Complex

The clique complex (flag complex) Cl(G) has:
- k-simplex σ = {v₀, ..., vₖ} iff the induced subgraph G[σ] is a complete graph K_{k+1}

This means every complete subgraph "fills in" to become a simplex. A triangle of edges becomes a filled triangle (2-simplex). A K₄ becomes a tetrahedron (3-simplex).

### Simplicial Homology

Given a simplicial complex K, define chain groups Cₖ = free abelian group on k-simplices, with boundary operator:

```
∂ₖ: Cₖ → C_{k−1}
∂ₖ([v₀, ..., vₖ]) = Σᵢ (−1)ⁱ [v₀, ..., v̂ᵢ, ..., vₖ]
```

The k-th homology group: Hₖ = ker(∂ₖ) / im(∂ₖ₊₁)

Betti numbers: βₖ = rank(Hₖ)

### Euler Characteristic

For a simplicial complex:

```
χ = Σₖ (−1)ᵏ fₖ = Σₖ (−1)ᵏ βₖ
```

Where fₖ = number of k-simplices. This is the Euler-Poincaré formula: alternating sum of face counts equals alternating sum of Betti numbers. It's a topological invariant — doesn't depend on the triangulation.

### Graph Laplacian

**Combinatorial Laplacian**: L = D − A where D is the degree matrix and A is the adjacency matrix.

```
L_{ij} = deg(i)     if i = j
L_{ij} = −1         if (i,j) is an edge
L_{ij} = 0          otherwise
```

**Normalized Laplacian**: L_norm = D^{−1/2} L D^{−1/2}

Key spectral properties:
- λ₁ = 0 always (constant vector is eigenvector)
- Multiplicity of λ = 0 equals β₀ (number of components)
- **Algebraic connectivity** λ₂ > 0 iff graph is connected
- λ₂ small → loosely connected, slow mixing
- λ₂ large → tightly connected, fast mixing

## Installation

```bash
cargo add graph-homology
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
graph-homology = "0.1"
```

## Related Crates

Part of the **SuperInstance Exocortex** math fleet:

- **[cohomology-ring](https://github.com/SuperInstance/cohomology-ring)** — Cup products and cohomology operations
- **[sheaf-laplacian](https://github.com/SuperInstance/sheaf-laplacian)** — Sheaf Laplacian and Hodge decomposition
- **[persistent-agent](https://github.com/SuperInstance/persistent-agent)** — Persistent homology for agent behavior
- **[tropical-graph](https://github.com/SuperInstance/tropical-graph)** — Max-plus algebra on graphs
- **[markov-blanket](https://github.com/SuperInstance/markov-blanket)** — Statistical boundary detection

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project.
