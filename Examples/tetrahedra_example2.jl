using NewRivaraProductions

m = NewRivaraProductions.simpleTetrahedralMesh()
NewRivaraProductions.write_vtk(m, "tetrahedra_example2_0.vtu")

NewRivaraProductions.get_edges(m.root)[3].x.MR = true # Refine the edge [1, 1, 0] - [1, 1, 1]
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "tetrahedra_example2_1.vtu")

NewRivaraProductions.get_edges(m.root.x.next)[5].x.MR = true # Refine the edge [1, 1, 0.5] - [1, 1, 1]
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "tetrahedra_example2_2.vtu")

NewRivaraProductions.get_edges(m.root.x.next)[5].x.MR = true # Refine the edge [1, 1, 0.5] - [1, 1, 1]
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "tetrahedra_example2_3.vtu")
