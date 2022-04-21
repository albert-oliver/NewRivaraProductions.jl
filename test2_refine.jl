using NewRivaraProductions

function mark_all_triangles_for_refinement(m)
    for t in NewRivaraProductions.collect_all_triangles(m)
        t.MR = true
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

m = my_setup(20)

@time NewRivaraProductions.refine!(m)

# using BenchmarkTools
# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60
# b = @benchmark NewRivaraProductions.refine!(m) setup=(m = my_setup(20))
# io = IOBuffer()
# show(io, "text/plain", b)
# s = String(take!(io))
# println(s)
