module TriangularMeshBenchmarks
    
using NewRivaraProductions
using BenchmarkTools

suite = BenchmarkGroup()

function mark_all_triangles_for_refinement(m)
    for t in NewRivaraProductions.collect_all_elements(m)
        t.x.MR = true
    end
end

function create_mesh_and_refine_nsteps(nsteps)
    m = NewRivaraProductions.simpleTriangularMesh()
    for _ in 1:nsteps
        mark_all_triangles_for_refinement(m)
        NewRivaraProductions.refine!(m)
    end
    return m
end

nsteps = 21
nseconds = 225

suite["all_refinements"] = @benchmarkable create_mesh_and_refine_nsteps(nsteps) seconds=nseconds
suite["last_refinement"] = @benchmarkable NewRivaraProductions.refine!(m) setup=(m = create_mesh_and_refine_nsteps(nsteps-1); mark_all_triangles_for_refinement(m)) seconds=nseconds

end

return TriangularMeshBenchmarks.suite
