using NewRivaraProductions
using BenchmarkTools

function mark_all_triangles_for_refinement(m)
    for i in m.triangles
        m.all_triangles[i].MR = true
    end
end

function my_setup(nsteps)
    m = NewRivaraProductions.simpleMesh()
    mark_all_triangles_for_refinement(m)
    for i in 1:nsteps
        NewRivaraProductions.refine!(m)
        mark_all_triangles_for_refinement(m)
    end
    return m
end

my_setup(1)

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 50
b = @benchmark NewRivaraProductions.refine!(m) setup=(m = my_setup(14))
io = IOBuffer()
show(io, "text/plain", b)
s = String(take!(io))
println(s)
