using NewRivaraProductions
m = NewRivaraProductions.simpleMesh()
NewRivaraProductions.write_vtk(m, "kk0.vtu")
m.triangles[1].edges[1].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "kk1.vtu")
m.triangles[1].edges[2].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "kk2.vtu")
m.triangles[5].edges[3].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "kk3.vtu")
m.triangles[11].edges[3].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "kk4.vtu")
