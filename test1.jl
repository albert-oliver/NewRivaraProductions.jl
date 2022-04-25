using NewRivaraProductions
m = NewRivaraProductions.simpleMesh()
NewRivaraProductions.write_vtk(m, "kk0.vtu")

triangles = NewRivaraProductions.collect_all_triangles(m)
triangles[1].x.edges[1].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "kk1.vtu")

triangles = NewRivaraProductions.collect_all_triangles(m)
triangles[3].x.edges[2].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "kk2.vtu")

triangles = NewRivaraProductions.collect_all_triangles(m)
triangles[7].x.edges[1].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "kk3.vtu")

triangles = NewRivaraProductions.collect_all_triangles(m)
triangles[9].x.edges[3].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "kk4.vtu")
