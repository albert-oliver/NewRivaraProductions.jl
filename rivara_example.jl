using NewRivaraProductions

coords = transpose([
    2.0 11.5 0.0;
    0.0 10.0 0.0;
    -1.0 -2.0 0.0;
    5.0 -7.0 0.0;
    10.0 0.0 0.0;
    10.1 6.0 0.0;
    5.0 10.5 0.0])

conec = transpose([
    1 2 7;
    2 3 7;
    3 4 7;
    4 5 7;
    5 6 7])

m = NewRivaraProductions.Mesh(coords, conec)
NewRivaraProductions.write_vtk(m, "rivara_example_0.vtu")

m.triangles[1].MR = true
NewRivaraProductions.refine!(m)
NewRivaraProductions.write_vtk(m, "rivara_example_1.vtu")
