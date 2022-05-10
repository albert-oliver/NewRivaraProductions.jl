module NewRivaraProductions

using LinearAlgebra
using StaticArrays
using WriteVTK

# Write your package code here.

abstract type AbstractMesh end
abstract type AbstractElement end

struct Node
    uvw::SVector{3,Float64}
    xyz::SVector{3,Float64}
end

mutable struct Edge
    nodes::SVector{2,Base.RefValue{Node}}
    nodes_id::SVector{2,Int}
    MR::Bool
    @atomic BR::Bool
    @atomic NA::Int
    sons::Union{SVector{2,Base.RefValue{Edge}}, Nothing}
end

mutable struct Triangle <: AbstractElement
    edges::SVector{3,Base.RefValue{Edge}}
    MR::Bool
    prev::Union{Base.RefValue{Triangle}, Nothing}
    next::Union{Base.RefValue{Triangle}, Nothing}
    sons::Union{SVector{2,Base.RefValue{Triangle}}, Nothing}
end

mutable struct Tetrahedron <: AbstractElement
    faces::SVector{4,Base.RefValue{Triangle}}
    MR::Bool
    BR::Bool
    prev::Union{Base.RefValue{Tetrahedron}, Nothing}
    next::Union{Base.RefValue{Tetrahedron}, Nothing}
end

mutable struct TriangularMesh <: AbstractMesh
    root::Base.RefValue{Triangle}
    nodes::Vector{Base.RefValue{Node}}
end

mutable struct TetrahedralMesh <: AbstractMesh
    root::Base.RefValue{Tetrahedron}
    nodes::Vector{Base.RefValue{Node}}
end

function get_triangle_edges!(conec, edges_per_node, nodes)
    edges = Vector{Base.RefValue{Edge}}(undef, 3)

    for j in 1:3
        n1, n2 = conec[circshift(1:3, -j)[1:2]]
        edge = haskey(edges_per_node, n1) && haskey(edges_per_node, n2) ? intersect(edges_per_node[n1], edges_per_node[n2]) : []
        if isempty(edge)
            edge = Ref(Edge([nodes[n1], nodes[n2]], [n1, n2], false, false, 1, nothing))
            push!(get!(edges_per_node, n1, []), edge)
            push!(get!(edges_per_node, n2, []), edge)
        else
            edge = edge[1]
            @atomic edge.x.NA += 1
        end
        edges[j] = edge
    end

    return edges
end

function get_tetrahedron_triangles!(conec, triangles_per_node, edges_per_node, nodes)
    triangles = Vector{Base.RefValue{Triangle}}()
    sizehint!(triangles, 4)

    tri_conecs = transpose([2 4 3;
                            1 3 4;
                            2 1 4;
                            1 2 3])

    for tri_conec in eachcol(conec[tri_conecs])
        n1, n2, n3 = tri_conec
        triangle = haskey(triangles_per_node, n1) && haskey(triangles_per_node, n2) && haskey(triangles_per_node, n3) ?
            intersect(triangles_per_node[n1], triangles_per_node[n2], triangles_per_node[n3]) :
            []
        if isempty(triangle)
            edges = get_triangle_edges!(tri_conec[:], edges_per_node, nodes)
            triangle = Ref(Triangle(edges, false, nothing, nothing, nothing))
            push!(get!(triangles_per_node, n1, []), triangle)
            push!(get!(triangles_per_node, n2, []), triangle)
            push!(get!(triangles_per_node, n3, []), triangle)
        else
            triangle = triangle[1]
        end
        push!(triangles, triangle)
    end

        return triangles
end

function TriangularMesh(coords::AbstractArray{Float64}, conec::AbstractMatrix{Int})
    dim, nnodes = size(coords)
    nconec, nelem = size(conec)

    dim ≠ 3 && throw(ArgumentError("coords should be a matrix of size 3×nNodes"))
    nconec ≠ 3 && throw(ArgumentError("conec should be a matrix of size 3×nTriangles"))

    nodes = Vector{Base.RefValue{Node}}(undef, nnodes)

    for i in 1:nnodes
        nodes[i] = Ref(Node(coords[:,i], coords[:,i]))
    end

    triangles = Vector{Triangle}(undef, nelem)
    edges_per_node = Dict{Int, Vector{Base.RefValue{Edge}}}()

    for i in 1:nelem
        edges = get_triangle_edges!(conec[:,i], edges_per_node, nodes)
        triangles[i] = Triangle(edges, false, nothing, nothing, nothing)
    end

    for i in 1:Base.length(triangles)-1
        triangles[i].next = Ref(triangles[i+1])
    end

    for i in 2:Base.length(triangles)
        triangles[i].prev = Ref(triangles[i-1])
    end

    return TriangularMesh(Ref(triangles[1]), nodes)
end

function TetrahedralMesh(coords::AbstractArray{Float64}, conec::AbstractMatrix{Int})
    dim, nnodes = size(coords)
    nconec, nelem = size(conec)

    dim ≠ 3 && throw(ArgumentError("coords should be a matrix of size 3×nNodes"))
    nconec ≠ 4 && throw(ArgumentError("conec should be a matrix of size 4×nTetrahedra"))

    nodes = Vector{Base.RefValue{Node}}(undef, nnodes)

    for i in 1:nnodes
        nodes[i] = Ref(Node(coords[:,i], coords[:,i]))
    end

    tetrahedra = Vector{Tetrahedron}(undef, nelem)
    edges_per_node = Dict{Int, Vector{Base.RefValue{Edge}}}()
    triangles_per_node = Dict{Int, Vector{Base.RefValue{Triangle}}}()

    for i in 1:nelem
        triangles = get_tetrahedron_triangles!(conec[:,i], triangles_per_node, edges_per_node, nodes)
        tetrahedra[i] = Tetrahedron(triangles, false, false, nothing, nothing)
    end

    for i in 1:Base.length(tetrahedra)-1
        tetrahedra[i].next = Ref(tetrahedra[i+1])
    end

    for i in 2:Base.length(tetrahedra)
        tetrahedra[i].prev = Ref(tetrahedra[i-1])
    end

    return TetrahedralMesh(Ref(tetrahedra[1]), nodes)
end

function simpleTriangularMesh()
    coords = transpose([0.0 0.0 0.0;
                        1.0 0.0 0.0;
                        1.0 1.0 0.0;
                        0.0 1.0 0.0])

    conec = transpose([2 3 1
                       4 1 3])

    return TriangularMesh(coords, conec)
end

function simpleTetrahedralMesh()


    coords = transpose([0.0 0.0 0.0;
                        1.0 0.0 0.0;
                        1.0 1.0 0.0;
                        0.0 1.0 0.0;
                        0.0 0.0 1.0;
                        1.0 0.0 1.0;
                        1.0 1.0 1.0;
                        0.0 1.0 1.0])

    conec = transpose([2 3 1 7;
                       6 2 1 7;
                       5 6 1 7;
                       8 5 1 7;
                       4 8 1 7;
                       3 4 1 7])

    return TetrahedralMesh(coords, conec)
end

length(e::Base.RefValue{Edge}) = norm(e.x.nodes[1].x.xyz - e.x.nodes[2].x.xyz)
new_coords(e::Base.RefValue{Edge}) = ((e.x.nodes[1]).x.xyz + (e.x.nodes[2]).x.xyz)/2

function common_node(e1::Base.RefValue{Edge}, e2::Base.RefValue{Edge})
    for n1 in e1.x.nodes, n2 in e2.x.nodes
        if n1 == n2
            return [n1]
        end
    end
    return Vector{Base.RefValue{Node}}()
end

function common_node_id(e1::Base.RefValue{Edge}, e2::Base.RefValue{Edge})
    for n1 in e1.x.nodes_id, n2 in e2.x.nodes_id
        if n1 == n2
            return n1
        end
    end
    @assert false "common_node_id needs two edges that actually share a node"
end

function common_edges(tri1::Triangle, tri2::Triangle)
    for e1 in tri1.edges, e2 in tri2.edges
        if e1 == e2
            return [e1]
        end
    end
    return Vector{Base.RefValue{Edge}}()
end
common_edges(tri1::Base.RefValue{Triangle}, tri2::Base.RefValue{Triangle}) = common_edges(tri1.x, tri2.x)

have_one_common_edge(tri1::Triangle, tri2::Triangle) = Base.length(common_edges(tri1, tri2)) == 1
have_one_common_edge(tri1::Base.RefValue{Triangle}, tri2::Base.RefValue{Triangle}) = have_one_common_edge(tri1.x, tri2.x)

isbroken(tri::Triangle) = !isnothing(tri.sons)
isbroken(tri::Base.RefValue{Triangle}) = isbroken(tri.x)

isbroken(e::Edge) = @atomic e.BR
isbroken(e::Base.RefValue{Edge}) = isbroken(e.x)

isbroken(tet::Tetrahedron) = tet.BR
isbroken(tet::Base.RefValue{Tetrahedron}) = isbroken(tet.x)

canbebroken(e::Edge) = (@atomicreplace e.BR false => true)[2]
canbebroken(e::Base.RefValue{Edge}) = canbebroken(e.x)

# A triangle or a tetrahedron is nonconformal if any of its edges is broken
isnonconformal(tri::Triangle) = any(isbroken, get_edges(tri))
isnonconformal(tri::Base.RefValue{Triangle}) = isnonconformal(tri.x)

isnonconformal(tet::Tetrahedron) = any(isbroken, get_edges(tet))
isnonconformal(tet::Base.RefValue{Tetrahedron}) = isnonconformal(tet.x)

ismarkedforrefinement(e::Edge) = e.MR
ismarkedforrefinement(e::Base.RefValue{Edge}) = ismarkedforrefinement(e.x)

# Usually, to see if a triangle is marked to be refined it should check itself and its edges.
# But, since the algorithm first bisects all the edges that are marked to be refined, we can ignore this second condition since it will always be false
ismarkedforrefinement(tri::Triangle) = tri.MR # || any(ismarkedforrefinement, get_edges(tri))
ismarkedforrefinement(tri::Base.RefValue{Triangle}) = ismarkedforrefinement(tri.x)

ismarkedforrefinement(tet::Tetrahedron) = tet.MR || any(ismarkedforrefinement, get_triangles(tet))
ismarkedforrefinement(tet::Base.RefValue{Tetrahedron}) = ismarkedforrefinement(tet.x)

get_root(m::AbstractMesh) = m.root

function Base.isless(e1::Base.RefValue{Edge}, e2::Base.RefValue{Edge})

    # First we just compare lengths
    l1 = length(e1)
    l2 = length(e2)
    if l1 ≉ l2
        return l1 < l2
    end

    # Then, if any edge is marked to be refined (or already broken)
    br1 = isbroken(e1) || e1.x.MR
    br2 = isbroken(e1) || e2.x.MR
    if br1 ≠ br2
        # This function returns the edge that is NOT going to be refined,
        # we want to refine the edge that is marked for refinement or already broken,
        # therefore if it's marked or already broken it should return false
        return !(br1)
    end

    # Next, the number of adjacent elements
    na1 = e1.x.NA
    na2 = e2.x.NA
    if na1 ≠ na2
        # This function returns the edge that is NOT going to be refined,
        # we want to refine the edge that has less adjacent elements,
        # therefore, we return the edge that has more adjacent elements.
        return na1 > na2
    end

    # Finally, to make the comparison deterministic, we compare the positions
    # of the center of the edge (the new vertex)
    x1, y1, z1 = new_coords(e1)
    x2, y2, z2 = new_coords(e2)
    if x1 ≉ x2
        return x1 < x2
    elseif y1 ≉ y2
        return y1 < y2
    elseif z1 ≉ z2
        return z1 < z2
    end

    # The only situation that I can think of is when e1 and e2 are the same edge...
    @assert e1 == e2 "Something weird happened"

    # Therefore, they are equal so e1 is NOT less than e2
    return false
end

get_edges(tri::Triangle) = tri.edges
get_edges(tri::Base.RefValue{Triangle}) = get_edges(tri.x)

get_edges(tet::Base.RefValue{Tetrahedron}) = get_edges(tet.x)
function get_edges(tet::Tetrahedron)
    edges = Vector{Base.RefValue{Edge}}()

    for e in [tet.faces[1].x.edges tet.faces[2].x.edges tet.faces[3].x.edges]
        if !any(isequal(e), edges)
            push!(edges, e)
        end
    end
    return SVector{6, Base.RefValue{Edge}}(edges)
end

get_max_edge(element::AbstractElement) = maximum(get_edges(element))

get_triangles(tet::Tetrahedron) = tet.faces

function get_sorted_edges(triangle)
    edges = Vector{Base.RefValue{Edge}}(get_edges(triangle))
    e_idx = findmax(edges)[2]

    return circshift(edges, -(e_idx-1))
end

function get_sorted_faces(tetrahedron::Tetrahedron)
    e1 = get_max_edge(tetrahedron)

    face_contains_e1 = [any(isequal(e1), get_edges(tri)) for tri in get_triangles(tetrahedron)]

    if (face_contains_e1[1] && face_contains_e1[2])
        return tetrahedron.faces[[1, 2, 3, 4]]
    elseif (face_contains_e1[1] && face_contains_e1[3])
        return tetrahedron.faces[[3, 1, 2, 4]]
    elseif (face_contains_e1[1] && face_contains_e1[4])
        return tetrahedron.faces[[1, 4, 2, 3]]
    elseif (face_contains_e1[2] && face_contains_e1[3])
        return tetrahedron.faces[[2, 3, 1, 4]]
    elseif (face_contains_e1[2] && face_contains_e1[4])
        return tetrahedron.faces[[4, 2, 3, 1]]
    elseif (face_contains_e1[3] && face_contains_e1[4])
        return tetrahedron.faces[[3, 4, 1, 2]]
    end
end

lk = ReentrantLock()

function add_new_elements!(m::AbstractMesh, tp::T, rt1::Base.RefValue{T}, rt2::Base.RefValue{T}) where T <: AbstractElement
    lock(lk) do
        rt1.x.prev = tp.prev
        if !isnothing(rt1.x.prev)
            rt1.x.prev.x.next = rt1
        else
            m.root = rt1
        end
        rt1.x.next = rt2
        rt2.x.prev = rt1
        rt2.x.next = tp.next
        if !isnothing(rt2.x.next)
            rt2.x.next.x.prev = rt2
        end

        # Clean (probably not needed)
        tp.prev = nothing
        tp.next = nothing
    end
    return nothing
end

function prod_bisect_edges!(m::AbstractMesh, element::AbstractElement)

    # First of all we can bisect all the edges that are marked to be refined...
    any_bisection = false
    for e in get_edges(element)
        if e.x.MR && canbebroken(e)
            bisect_edge!(m, e)
            any_bisection = true
        end
    end

    # Then, we'll see if the longest edge needs to be refined (because the triangle is marked to be refined, or the element needs to be conformed)
    e1 = get_max_edge(element)

    if (ismarkedforrefinement(element) || isnonconformal(element)) && canbebroken(e1)
        bisect_edge!(m, e1)
        return true
    end

    return any_bisection

end

function bisect_edge!(m::AbstractMesh, edge::Base.RefValue{Edge})

    edge.x.MR = false

    # Generate new node
    new_node = Ref(Node(new_coords(edge), new_coords(edge)))

    new_node_id = Vector{Int}(undef, Threads.nthreads())

    # Add the new node in the common vector and getting the id.
    # We need to lock
    @lock lk (push!(m.nodes, new_node); new_node_id[Threads.threadid()] = Base.length(m.nodes))

    # Generate the first new edge
    edge1 = Ref(Edge([edge.x.nodes[1], new_node], [edge.x.nodes_id[1], new_node_id[Threads.threadid()]], false, false, edge.x.NA, nothing))

    # Generate the second new edge
    edge2 = Ref(Edge([new_node, edge.x.nodes[2]], [new_node_id[Threads.threadid()], edge.x.nodes_id[2]], false, false, edge.x.NA, nothing))

    # Add the sons of the initial edge
    edge.x.sons = [edge1, edge2]

    return nothing

end

function bisect_triangle!(triangle::Triangle)

    edge1, edge2, edge3 = get_sorted_edges(triangle)

    triangle.MR = false

    edge4, edge5 = edge1.x.sons

    if (isempty(common_node(edge3, edge4)))
        edge4, edge5 = edge5, edge4
    end

    v1 = common_node(edge4, edge5)[]
    v2 = common_node(edge2, edge3)[]
    v1_id = common_node_id(edge4, edge5)
    v2_id = common_node_id(edge2, edge3)
    new_edge = Ref(Edge([v1, v2], [v1_id, v2_id], false, false, 2, nothing))

    triangle1 = Ref(Triangle([edge3, edge4, new_edge], false, nothing, nothing, nothing))
    triangle2 = Ref(Triangle([edge2, new_edge, edge5], false, nothing, nothing, nothing))

    triangle.sons = [triangle1, triangle2]

    return triangle1, triangle2
end

function bisect_tetrahedron!(tetrahedron::Tetrahedron, sorted_faces)

    tetrahedron.MR = false

    face1, face2, face3, face4 = sorted_faces

    face5, face7 = face1.x.sons
    if !have_one_common_edge(face3, face5)
        face5, face7 = face7, face5
    end

    face6, face8 = face2.x.sons
    if !have_one_common_edge(face3, face6)
        face6, face8 = face8, face6
    end

    # New face
    e1 = common_edges(face3, face4)[]
    e2 = common_edges(face6, face8)[]
    e3 = common_edges(face5, face7)[]

    @atomic e1.x.NA += 1
    @atomic e2.x.NA += 1
    @atomic e3.x.NA += 1

    new_face = Ref(Triangle([e1, e2, e3], false, nothing, nothing, nothing))

    tet1 = Ref(Tetrahedron([face4, face7, new_face, face8], false, false, nothing, nothing))
    tet2 = Ref(Tetrahedron([face3, face6, new_face, face5], false, false, nothing, nothing))

    return tet1, tet2

end

function prod_bisect_element!(m::AbstractMesh, triangle::Triangle)

    if !isbroken(get_max_edge(triangle))
        return false
    end

    rt1, rt2 = bisect_triangle!(triangle)

    add_new_elements!(m, triangle, rt1, rt2)

    prod_bisect_element!(m, rt1.x)
    prod_bisect_element!(m, rt2.x)

    return true
end

function prod_bisect_element!(m::AbstractMesh, tetrahedron::Tetrahedron)

    if !isbroken(get_max_edge(tetrahedron))
        return false
    end

    sorted_faces = get_sorted_faces(tetrahedron)

    if !isbroken(sorted_faces[1])
        lock(lk) do # We don't want the other adjacent tetrahedron to break the triangle at the same time
            bisect_triangle!(sorted_faces[1].x)
        end
    end

    if !isbroken(sorted_faces[2])
        lock(lk) do # We don't want the other adjacent tetrahedron to break the triangle at the same time
            bisect_triangle!(sorted_faces[2].x)
        end
    end
    
    tet1, tet2 = bisect_tetrahedron!(tetrahedron, sorted_faces)

    add_new_elements!(m, tetrahedron, tet1, tet2)

    prod_bisect_element!(m, tet1.x)
    prod_bisect_element!(m, tet2.x)

    return true

end

function collect_all_elements(m::AbstractMesh)
    node = get_root(m)
    elements = Vector{typeof(node)}()
    while !isnothing(node)
        push!(elements, node)
        node = node.x.next
    end
    return elements
end

function refine!(m::AbstractMesh)

    run = true

    while run
        run = false

        l_triangles = collect_all_elements(m)

        edges_broken = trues(Threads.nthreads())
        while any(edges_broken)
            edges_broken .= false
            Threads.@threads for t in l_triangles
                edges_broken[Threads.threadid()] |= prod_bisect_edges!(m, t.x)
            end
            run |= any(edges_broken)
        end

        if run
            Threads.@threads for t in l_triangles
                prod_bisect_element!(m, t.x)
            end
        end
    end

    return nothing

end

collect_all_nodes(m::AbstractMesh) = m.nodes
get_conec(t::Base.RefValue{Triangle}) = get_conec(t.x)
get_conec(t::Triangle) = [
    common_node_id(t.edges[2], t.edges[3]),
    common_node_id(t.edges[3], t.edges[1]),
    common_node_id(t.edges[1], t.edges[2])
]

function three_vector_intersect(a, b, c)
    for m in a, n in b, t in c
        if m == n == t
            return m
        end
    end
    @assert false "three_vector_intersect needs three vectors that actually share an element"
end

function get_conec(t::Tetrahedron)
    nodes_ids = SMatrix{3,4,Int}(vcat([get_conec(f) for f in t.faces]...))

    return SVector{4,Int}(
        three_vector_intersect(nodes_ids[:,2], nodes_ids[:,3], nodes_ids[:,4]),
        three_vector_intersect(nodes_ids[:,3], nodes_ids[:,4], nodes_ids[:,1]),
        three_vector_intersect(nodes_ids[:,4], nodes_ids[:,1], nodes_ids[:,2]),
        three_vector_intersect(nodes_ids[:,1], nodes_ids[:,2], nodes_ids[:,3])
    )
end

get_VTKCellType(_::TriangularMesh) = VTKCellTypes.VTK_TRIANGLE
get_VTKCellType(_::TetrahedralMesh) = VTKCellTypes.VTK_TETRA

function get_xyz_uvw(m::AbstractMesh)
    nodes_vec = collect_all_nodes(m)
    xyz = similar(nodes_vec[1].x.xyz, 3, Base.length(nodes_vec))
    uvw = similar(xyz)
    Threads.@threads for i in 1:Base.length(nodes_vec)
        xyz[:,i] .= nodes_vec[i].x.xyz
        uvw[:,i] .= nodes_vec[i].x.uvw
    end
    return xyz, uvw
    # xyz = mapreduce(n -> n.x.xyz, hcat, nodes_vec)
    # uvw = mapreduce(n -> n.x.uvw, hcat, nodes_vec)
    # return xyz, uvw
end

function get_VTKconecs(m::AbstractMesh) 
    vtk_cell_type = get_VTKCellType(m)
    tetra = collect_all_elements(m)
    vtkc1 = [MeshCell(vtk_cell_type, get_conec(tetra[1].x))]
    vtk_conecs = similar(vtkc1,  Base.length(tetra))
    Threads.@threads for i in 1:Base.length(vtk_conecs)
        vtk_conecs[i] = MeshCell(vtk_cell_type, get_conec(tetra[i].x))
    end
    return vtk_conecs
    # return [MeshCell(vtk_cell_type, get_conec(t.x)) for t in collect_all_elements(m)]
end

function write_vtk(m::AbstractMesh, filename)
    coords, uvw = get_xyz_uvw(m)
    vtk_conecs = get_VTKconecs(m)
    vtk_grid(filename, coords, vtk_conecs) do vtk
        vtk["uvw"] = uvw
    end
end

end
