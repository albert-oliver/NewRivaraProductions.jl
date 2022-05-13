using NewRivaraProductions
m = NewRivaraProductions.simpleTriangularMesh()
NewRivaraProductions.write_vtk(m, "triangle_example_0.vtu")

triangles = NewRivaraProductions.collect_all_elements(m)
triangles[1].edges[2].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "triangle_example_1.vtu")

triangles = NewRivaraProductions.collect_all_elements(m)
triangles[1].edges[2].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "triangle_example_2.vtu")

triangles = NewRivaraProductions.collect_all_elements(m)
triangles[7].edges[1].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "triangle_example_3.vtu")

triangles = NewRivaraProductions.collect_all_elements(m)
triangles[13].edges[3].x.MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "triangle_example_4.vtu")
