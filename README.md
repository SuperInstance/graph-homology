# graph-homology

> **The topology hidden in graphs. Clique complexes, graph Laplacians, Euler characteristic.**

[![crates.io](https://img.shields.io/crates/v/graph-homology.svg)](https://crates.io/crates/graph-homology)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Computes homology of graph networks via clique complexes. Every complete subgraph becomes a simplex, giving the graph a topological structure. From this, we extract Betti numbers, Euler characteristic, and graph Laplacians.

## The Clique Complex

Given a graph G, the **clique complex** Cl(G) is the simplicial complex whose simplices are the complete subgraphs of G:
- Nodes → 0-simplices (vertices)
- Edges → 1-simplices
- Triangles → 2-simplices
- K₄ subgraphs → 3-simplices
- ...

This turns graph theory into topology.

## What It Computes

- **H₀**: connected components
- **H₁**: independent cycles (graph cycles not filled by triangles)
- **Euler characteristic**: V - E + F
- **Combinatorial Laplacian**: L = D - A
- **Normalized Laplacian**: D^{-1/2} L D^{-1/2}

## License

MIT © [SuperInstance](https://github.com/SuperInstance)
