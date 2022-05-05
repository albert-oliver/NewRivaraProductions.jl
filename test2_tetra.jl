using NewRivaraProductions

m = NewRivaraProductions.simpleTetrahedralMesh()
NewRivaraProductions.write_vtk(m, "test2_tetra_0.vtu")

NewRivaraProductions.get_edges(m.root)[1].x.MR = true # Refine the edge [1, 1, 0] - [1, 1, 1]
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "test2_tetra_1.vtu")

NewRivaraProductions.get_edges(m.root)[2].x.MR = true # Refine the edge [1, 1, 0.5] - [1, 1, 1]
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "test2_tetra_2.vtu")

NewRivaraProductions.get_edges(m.root)[2].x.MR = true # Refine the edge [1, 1, 0.75] - [1, 1, 1]
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "test2_tetra_3.vtu")
