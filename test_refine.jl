using NewRivaraProductions

function mark_all_triangles_for_refinement(m)
    for t in NewRivaraProductions.collect_all_triangles(m)
        t.x.MR = true
    end
end

function test_refine(nsteps)
    m = NewRivaraProductions.simpleMesh()
    for i in 1:nsteps
        mark_all_triangles_for_refinement(m)
        NewRivaraProductions.refine!(m)
    end
end

test_refine(1)

@time test_refine(21)
