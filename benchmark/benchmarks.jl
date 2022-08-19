using BenchmarkTools

const SUITE = BenchmarkGroup()

SUITE["TetrahedralMesh"] = include("./benchmark_tetrahedra.jl")
SUITE["TriangularMesh"] = include("./benchmark_triangles.jl")
