using NewRivaraProductions

function mark_all_triangles_for_refinement(m)
    for i in m.triangles
        m.all_triangles[i].MR = true
    end
end

function test_refine()
    m = NewRivaraProductions.simpleMesh()
    for i in 1:15
        mark_all_triangles_for_refinement(m)
        NewRivaraProductions.refine!(m)
    end
end

test_refine()

@time test_refine()

