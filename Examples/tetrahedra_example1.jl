using NewRivaraProductions

coords = transpose([0.0 0.0 0.0
                    1.0 1.0 1.0
                    1.0 1.0 0.0
                    1.0 0.0 0.0])

conec = transpose([2 1 3 4])

m = NewRivaraProductions.TetrahedralMesh(coords, conec)
NewRivaraProductions.write_vtk(m, "tetrahedra_example1_0.vtu")

m.root.x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "tetrahedra_example1_1.vtu")

m.root.x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "tetrahedra_example1_2.vtu")

m.root.x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "tetrahedra_example1_3.vtu")
