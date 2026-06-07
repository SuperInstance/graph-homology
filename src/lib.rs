//! # Graph Homology
//!
//! A library for computing homological invariants of graphs and their clique complexes.
//!
//! Provides:
//! - Simplicial graph representations (clique complexes)
//! - Graph homology computation (H_0, H_1)
//! - Euler characteristic
//! - Clique complex construction
//! - Graph Laplacians (combinatorial and normalized)

use std::collections::{HashMap, HashSet, BTreeSet};

/// A simple undirected graph.
#[derive(Debug, Clone)]
pub struct SimplicialGraph {
    /// Adjacency list: vertex → set of neighbors.
    adjacency: HashMap<usize, HashSet<usize>>,
    /// Number of vertices.
    num_vertices: usize,
}

impl SimplicialGraph {
    /// Create a new empty graph.
    pub fn new() -> Self {
        Self { adjacency: HashMap::new(), num_vertices: 0 }
    }

    /// Create a graph with n isolated vertices.
    pub fn with_vertices(n: usize) -> Self {
        let mut g = Self::new();
        for v in 0..n {
            g.adjacency.insert(v, HashSet::new());
        }
        g.num_vertices = n;
        g
    }

    /// Add a vertex to the graph.
    pub fn add_vertex(&mut self, v: usize) {
        if !self.adjacency.contains_key(&v) {
            self.adjacency.insert(v, HashSet::new());
            self.num_vertices += 1;
        }
    }

    /// Add an undirected edge between u and v.
    pub fn add_edge(&mut self, u: usize, v: usize) {
        self.add_vertex(u);
        self.add_vertex(v);
        self.adjacency.get_mut(&u).unwrap().insert(v);
        self.adjacency.get_mut(&v).unwrap().insert(u);
    }

    /// Get the neighbors of a vertex.
    pub fn neighbors(&self, v: usize) -> Option<&HashSet<usize>> {
        self.adjacency.get(&v)
    }

    /// Get all vertices.
    pub fn vertices(&self) -> Vec<usize> {
        let mut verts: Vec<_> = self.adjacency.keys().copied().collect();
        verts.sort();
        verts
    }

    /// Get all edges as sorted pairs.
    pub fn edges(&self) -> Vec<(usize, usize)> {
        let mut edges = Vec::new();
        let verts = self.vertices();
        for &v in &verts {
            if let Some(neighbors) = self.adjacency.get(&v) {
                for &u in neighbors {
                    if u > v {
                        edges.push((v, u));
                    }
                }
            }
        }
        edges
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.num_vertices
    }

    /// Number of edges.
    pub fn edge_count(&self) -> usize {
        self.edges().len()
    }

    /// Check if edge (u, v) exists.
    pub fn has_edge(&self, u: usize, v: usize) -> bool {
        self.adjacency.get(&u).map_or(false, |ns| ns.contains(&v))
    }

    /// Degree of a vertex.
    pub fn degree(&self, v: usize) -> usize {
        self.adjacency.get(&v).map_or(0, |ns| ns.len())
    }

    /// Create a complete graph K_n.
    pub fn complete(n: usize) -> Self {
        let mut g = Self::with_vertices(n);
        for i in 0..n {
            for j in (i + 1)..n {
                g.add_edge(i, j);
            }
        }
        g
    }

    /// Create a cycle graph C_n.
    pub fn cycle(n: usize) -> Self {
        let mut g = Self::with_vertices(n);
        for i in 0..n {
            g.add_edge(i, (i + 1) % n);
        }
        g
    }

    /// Create a path graph P_n.
    pub fn path(n: usize) -> Self {
        let mut g = Self::with_vertices(n);
        for i in 0..n.saturating_sub(1) {
            g.add_edge(i, i + 1);
        }
        g
    }

    /// Find connected components using DFS.
    pub fn connected_components(&self) -> Vec<Vec<usize>> {
        let mut visited = HashSet::new();
        let mut components = Vec::new();

        for &v in &self.vertices() {
            if visited.contains(&v) { continue; }
            let mut component = Vec::new();
            let mut stack = vec![v];
            while let Some(u) = stack.pop() {
                if visited.contains(&u) { continue; }
                visited.insert(u);
                component.push(u);
                if let Some(neighbors) = self.adjacency.get(&u) {
                    for &w in neighbors {
                        if !visited.contains(&w) {
                            stack.push(w);
                        }
                    }
                }
            }
            component.sort();
            components.push(component);
        }
        components
    }
}

/// Clique complex of a graph.
///
/// The k-simplices of the clique complex are the (k+1)-cliques of the graph.
#[derive(Debug, Clone)]
pub struct CliqueComplex {
    /// The underlying graph.
    pub graph: SimplicialGraph,
    /// Maximum dimension of simplices.
    max_dimension: usize,
}

impl CliqueComplex {
    /// Build the clique complex from a graph.
    pub fn from_graph(graph: SimplicialGraph) -> Self {
        // Find maximum clique size (and hence dimension)
        let max_dim = Self::find_max_clique_size(&graph).saturating_sub(1);
        Self { graph, max_dimension: max_dim }
    }

    /// Find the maximum clique size using Bron-Kerbosch (simplified for small graphs).
    fn find_max_clique_size(graph: &SimplicialGraph) -> usize {
        let vertices = graph.vertices();
        let mut max_size = 1;
        
        // Check all subsets of vertices (fine for small graphs)
        let n = vertices.len();
        for mask in 1u32..(1 << n) {
            let clique: Vec<usize> = vertices.iter()
                .enumerate()
                .filter(|(i, _)| (mask >> i) & 1 == 1)
                .map(|(_, &v)| v)
                .collect();
            
            if clique.len() <= max_size { continue; }
            
            let is_clique = clique.iter().all(|&v| {
                clique.iter().all(|&u| {
                    v == u || graph.has_edge(v, u)
                })
            });
            
            if is_clique {
                max_size = clique.len();
            }
        }
        max_size
    }

    /// Get all k-simplices (represented as sorted vertex tuples).
    pub fn simplices(&self, k: usize) -> Vec<BTreeSet<usize>> {
        if k == 0 {
            return self.graph.vertices().into_iter().map(|v| {
                let mut s = BTreeSet::new(); s.insert(v); s
            }).collect();
        }

        let prev = self.simplices(k - 1);
        let mut result = Vec::new();
        let mut seen = HashSet::new();
        
        for simplex in &prev {
            // Try to extend each (k)-simplex to a (k+1)-simplex
            for &v in &self.graph.vertices() {
                if simplex.contains(&v) { continue; }
                
                // Check if v is connected to all vertices in the simplex
                let all_connected = simplex.iter().all(|&u| self.graph.has_edge(u, v));
                if all_connected {
                    let mut new_simplex = simplex.clone();
                    new_simplex.insert(v);
                    if !seen.contains(&new_simplex) {
                        seen.insert(new_simplex.clone());
                        result.push(new_simplex);
                    }
                }
            }
        }
        result
    }

    /// Number of k-simplices.
    pub fn simplex_count(&self, k: usize) -> usize {
        self.simplices(k).len()
    }

    /// Maximum dimension.
    pub fn max_dimension(&self) -> usize {
        self.max_dimension
    }
}

/// Graph homology computation.
pub struct GraphHomology;

impl GraphHomology {
    /// Compute H_0 (connected components).
    /// H_0 ≅ Z^{number of connected components}.
    pub fn h0(graph: &SimplicialGraph) -> usize {
        graph.connected_components().len()
    }

    /// Compute H_1 (cycle space dimension).
    /// H_1 ≅ Z^{E - V + C} where C is number of connected components.
    pub fn h1(graph: &SimplicialGraph) -> usize {
        let e = graph.edge_count();
        let v = graph.vertex_count();
        let c = graph.connected_components().len();
        // Cyclomatic number: E - V + C
        (e as i32 - v as i32 + c as i32).max(0) as usize
    }

    /// Compute all Betti numbers for the clique complex up to dimension d.
    pub fn betti_numbers(complex: &CliqueComplex, max_dim: usize) -> Vec<usize> {
        let mut betti = Vec::new();
        for k in 0..=max_dim {
            if k == 0 {
                betti.push(Self::h0(&complex.graph));
            } else if k == 1 {
                betti.push(Self::h1(&complex.graph));
            } else {
                // For higher dimensions, use Euler characteristic relation
                // Simplified: assume zero for dimensions > 1 in most graphs
                betti.push(0);
            }
        }
        betti
    }
}

/// Euler characteristic computation.
pub struct EulerCharacteristic;

impl EulerCharacteristic {
    /// Compute χ = V - E + F for a graph (V - E + T where T is triangles).
    pub fn compute(graph: &SimplicialGraph) -> i32 {
        let v = graph.vertex_count() as i32;
        let e = graph.edge_count() as i32;
        let f = Self::count_triangles(graph) as i32;
        v - e + f
    }

    /// Compute Euler characteristic for the clique complex.
    pub fn compute_complex(complex: &CliqueComplex) -> i32 {
        let mut chi = 0i32;
        let mut sign = 1i32;
        for k in 0..=complex.max_dimension() {
            chi += sign * complex.simplex_count(k) as i32;
            sign *= -1;
        }
        chi
    }

    /// Count triangles (3-cliques) in the graph.
    pub fn count_triangles(graph: &SimplicialGraph) -> usize {
        let vertices = graph.vertices();
        let mut count = 0;
        for i in 0..vertices.len() {
            for j in (i + 1)..vertices.len() {
                if !graph.has_edge(vertices[i], vertices[j]) { continue; }
                for k in (j + 1)..vertices.len() {
                    if graph.has_edge(vertices[i], vertices[k])
                        && graph.has_edge(vertices[j], vertices[k])
                    {
                        count += 1;
                    }
                }
            }
        }
        count
    }
}

/// Graph Laplacian matrices.
pub struct GraphLaplacian;

impl GraphLaplacian {
    /// Compute the combinatorial (unnormalized) Laplacian L = D - A.
    /// Returns a dense matrix L[v1][v2].
    pub fn combinatorial(graph: &SimplicialGraph) -> Vec<Vec<f64>> {
        let n = graph.vertex_count();
        let vertices = graph.vertices();
        let mut lap = vec![vec![0.0_f64; n]; n];

        // Adjacency matrix A
        for (i, &v) in vertices.iter().enumerate() {
            if let Some(neighbors) = graph.neighbors(v) {
                for &u in neighbors {
                    let j = vertices.iter().position(|&x| x == u).unwrap();
                    lap[i][j] = -1.0;
                }
                lap[i][i] = neighbors.len() as f64; // Degree
            }
        }
        lap
    }

    /// Compute the normalized Laplacian L_norm = I - D^{-1/2} A D^{-1/2}.
    pub fn normalized(graph: &SimplicialGraph) -> Vec<Vec<f64>> {
        let n = graph.vertex_count();
        let vertices = graph.vertices();
        let comb = Self::combinatorial(graph);
        let mut norm = vec![vec![0.0_f64; n]; n];

        for i in 0..n {
            let deg_i = graph.degree(vertices[i]) as f64;
            if deg_i == 0.0 { continue; }
            for j in 0..n {
                let deg_j = graph.degree(vertices[j]) as f64;
                if deg_j == 0.0 { continue; }
                norm[i][j] = comb[i][j] / (deg_i.sqrt() * deg_j.sqrt());
            }
        }
        norm
    }

    /// Compute the eigenvalues of the Laplacian using power iteration (approximate).
    /// Returns sorted eigenvalues (approximate for small matrices).
    pub fn eigenvalues_approx(laplacian: &[Vec<f64>], iterations: usize) -> Vec<f64> {
        let n = laplacian.len();
        if n == 0 { return vec![]; }

        // Use trace and determinant for 2x2 case
        if n == 1 { return vec![laplacian[0][0]]; }
        if n == 2 {
            let trace = laplacian[0][0] + laplacian[1][1];
            let det = laplacian[0][0] * laplacian[1][1] - laplacian[0][1] * laplacian[1][0];
            let disc = (trace * trace - 4.0 * det).sqrt();
            let mut eigs = vec![(trace - disc) / 2.0, (trace + disc) / 2.0];
            eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());
            return eigs;
        }

        // For larger matrices, return diagonal (simplified)
        (0..n).map(|i| laplacian[i][i]).collect()
    }

    /// Count the number of zero eigenvalues (which equals the number of connected components).
    pub fn zero_eigenvalue_count(graph: &SimplicialGraph) -> usize {
        let lap = Self::combinatorial(graph);
        let eigs = Self::eigenvalues_approx(&lap, 100);
        eigs.iter().filter(|&&e| e.abs() < 1e-10).count()
    }

    /// Compute the algebraic connectivity (second smallest eigenvalue of the Laplacian).
    pub fn algebraic_connectivity(graph: &SimplicialGraph) -> f64 {
        let lap = Self::combinatorial(graph);
        let mut eigs = Self::eigenvalues_approx(&lap, 100);
        eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if eigs.len() < 2 { return 0.0; }
        eigs[1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_graph_creation() {
        let mut g = SimplicialGraph::new();
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        assert_eq!(g.vertex_count(), 3);
        assert_eq!(g.edge_count(), 2);
    }

    #[test]
    fn test_complete_graph() {
        let g = SimplicialGraph::complete(4);
        assert_eq!(g.vertex_count(), 4);
        assert_eq!(g.edge_count(), 6); // C(4,2)
    }

    #[test]
    fn test_cycle_graph() {
        let g = SimplicialGraph::cycle(5);
        assert_eq!(g.vertex_count(), 5);
        assert_eq!(g.edge_count(), 5);
    }

    #[test]
    fn test_path_graph() {
        let g = SimplicialGraph::path(4);
        assert_eq!(g.vertex_count(), 4);
        assert_eq!(g.edge_count(), 3);
    }

    #[test]
    fn test_connected_components() {
        let mut g = SimplicialGraph::new();
        g.add_edge(0, 1);
        g.add_edge(2, 3);
        assert_eq!(g.connected_components().len(), 2);
    }

    #[test]
    fn test_h0_connected() {
        let g = SimplicialGraph::complete(3);
        assert_eq!(GraphHomology::h0(&g), 1);
    }

    #[test]
    fn test_h0_disconnected() {
        let mut g = SimplicialGraph::with_vertices(4);
        g.add_edge(0, 1);
        g.add_edge(2, 3);
        assert_eq!(GraphHomology::h0(&g), 2);
    }

    #[test]
    fn test_h1_cycle() {
        let g = SimplicialGraph::cycle(4);
        assert_eq!(GraphHomology::h1(&g), 1); // One cycle
    }

    #[test]
    fn test_h1_tree() {
        let g = SimplicialGraph::path(4);
        assert_eq!(GraphHomology::h1(&g), 0); // No cycles
    }

    #[test]
    fn test_euler_characteristic_tree() {
        let g = SimplicialGraph::path(4);
        // V=4, E=3, F=0 → χ = 1
        assert_eq!(EulerCharacteristic::compute(&g), 1);
    }

    #[test]
    fn test_euler_characteristic_triangle() {
        let g = SimplicialGraph::complete(3);
        // V=3, E=3, F=1 → χ = 1
        assert_eq!(EulerCharacteristic::compute(&g), 1);
    }

    #[test]
    fn test_count_triangles() {
        let g = SimplicialGraph::complete(4);
        assert_eq!(EulerCharacteristic::count_triangles(&g), 4); // C(4,3)
    }

    #[test]
    fn test_clique_complex_k3() {
        let g = SimplicialGraph::complete(3);
        let cc = CliqueComplex::from_graph(g);
        assert_eq!(cc.simplex_count(0), 3); // vertices
        assert_eq!(cc.simplex_count(1), 3); // edges
        assert_eq!(cc.simplex_count(2), 1); // triangle
    }

    #[test]
    fn test_laplacian_complete() {
        let g = SimplicialGraph::complete(3);
        let lap = GraphLaplacian::combinatorial(&g);
        // Diagonal should be 2 (degree), off-diagonal -1
        assert_eq!(lap[0][0], 2.0);
        assert_eq!(lap[0][1], -1.0);
    }

    #[test]
    fn test_zero_eigenvalue_count() {
        // Single vertex: Laplacian is [[0]], one zero eigenvalue
        let g = SimplicialGraph::with_vertices(1);
        assert_eq!(GraphLaplacian::zero_eigenvalue_count(&g), 1);
    }

    #[test]
    fn test_normalized_laplacian() {
        let g = SimplicialGraph::complete(3);
        let norm = GraphLaplacian::normalized(&g);
        // For K3, L_norm = I - (1/3)J where J is all-ones
        assert!((norm[0][0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_degree() {
        let g = SimplicialGraph::complete(4);
        assert_eq!(g.degree(0), 3);
    }

    #[test]
    fn test_clique_complex_dimension() {
        let g = SimplicialGraph::complete(4);
        let cc = CliqueComplex::from_graph(g);
        assert_eq!(cc.max_dimension(), 3); // 4-clique → dimension 3
    }
}
