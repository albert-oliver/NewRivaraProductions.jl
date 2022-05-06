using NewRivaraProductions

println("Number of Threads: $(Threads.nthreads())")

function mark_all_tetrahedra_for_refinement(m)
    for t in NewRivaraProductions.collect_all_elements(m)
        t.x.MR = true
    end
end

function test_refine(nsteps)
    m = NewRivaraProductions.simpleTetrahedralMesh()
    for i in 1:nsteps
        mark_all_tetrahedra_for_refinement(m)
        NewRivaraProductions.refine!(m)
    end
    return m
end

test_refine(1)

@time test_refine(18)
