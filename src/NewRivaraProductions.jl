module NewRivaraProductions

using LinearAlgebra
using StaticArrays
using WriteVTK

# Write your package code here.

struct Node
    uvw::SVector{3,Float64}
    xyz::SVector{3,Float64}
end

mutable struct Edge
    nodes::SVector{2,Base.RefValue{Node}}
    nodes_id::SVector{2,Int}
    MR::Bool
    NA::Int
    sons::Union{SVector{2,Base.RefValue{Edge}}, Nothing}
end

mutable struct Triangle
    edges::SVector{3,Base.RefValue{Edge}}
    MR::Bool
    BR::Bool
    prev::Union{Base.RefValue{Triangle}, Nothing}
    next::Union{Base.RefValue{Triangle}, Nothing}
end

mutable struct Mesh
    root_triangle::Base.RefValue{Triangle}
    nodes::Vector{Base.RefValue{Node}}
end

function Mesh(coords::AbstractArray{Float64}, conec::AbstractMatrix{Int})
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
        edges = Vector{Base.RefValue{Edge}}(undef, 3)

        for j in 1:3
            n1 = conec[j, i]
            n2 = conec[(j%3)+1, i]
            edge = haskey(edges_per_node, n1) && haskey(edges_per_node, n2) ? intersect(edges_per_node[n1], edges_per_node[n2]) : []
            if isempty(edge)
                edge = Ref(Edge([nodes[n1], nodes[n2]], [n1, n2], false, 1, nothing))
                push!(get!(edges_per_node, n1, []), edge)
                push!(get!(edges_per_node, n2, []), edge)
            else
                edge = edge[1]
                edge.x.NA = edge.x.NA + 1
            end
            edges[j] = edge
        end
        triangles[i] = Triangle(edges, false, false, nothing, nothing)
    end

    for i in 1:Base.length(triangles)-1
        triangles[i].next = Ref(triangles[i+1])
    end

    for i in 2:Base.length(triangles)
        triangles[i].prev = Ref(triangles[i-1])
    end

    return Mesh(Ref(triangles[1]), nodes)
end

function simpleMesh()
    coords = transpose([0.0 0.0 0.0;
                        1.0 0.0 0.0;
                        1.0 1.0 0.0;
                        0.0 1.0 0.0])

    conec = transpose([1 2 3;
                       3 4 1])

    return Mesh(coords, conec)
end

length(e::Base.RefValue{Edge}) = norm(e.x.nodes[1].x.xyz - e.x.nodes[2].x.xyz)
new_coords(e::Base.RefValue{Edge}) = ((e.x.nodes[1]).x.xyz + (e.x.nodes[2]).x.xyz)/2

common_node(e1::Base.RefValue{Edge}, e2::Base.RefValue{Edge}) = intersect(e1.x.nodes, e2.x.nodes)
common_node_id(e1::Base.RefValue{Edge}, e2::Base.RefValue{Edge}) = intersect(e1.x.nodes_id, e2.x.nodes_id)

isbroken(t::Triangle) = t.BR
isbroken(e::Base.RefValue{Edge}) = !isnothing(e.x.sons)

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

function get_edges(t::Triangle)
    return t.edges
end

function get_sorted_edges(triangle)
    edges = Vector{Base.RefValue{Edge}}(get_edges(triangle))
    e_idx = findmax(edges)[2]

    return circshift(edges, -(e_idx-1))
end

lk = ReentrantLock()

function add_new_triangles!(m::Mesh, tp::Triangle, t1::Triangle, t2::Triangle)
    lock(lk) do
        rt1 = Ref(t1)
        rt2 = Ref(t2)
        t1.prev = tp.prev
        if !isnothing(t1.prev)
            t1.prev.x.next = rt1
        else
            m.root_triangle = rt1
        end
        t1.next = rt2
        t2.prev = rt1
        t2.next = tp.next
        if !isnothing(t2.next)
            t2.next.x.prev = rt2
        end

        # Clean (probably not needed)
        tp.prev = nothing
        tp.next = nothing
    end
end

function p1_mark_edges!(m::Mesh, triangle::Base.RefValue{Triangle})

    if isbroken(triangle.x)
        return false
    end

    e1, e2, e3 = get_sorted_edges(triangle.x)

    if (triangle.x.MR || e1.x.MR || (isbroken(e2) || e2.x.MR) || (isbroken(e3) || e3.x.MR))
        lock(lk) do
            if (!isbroken(e1)) 
                p2_bisect_edge!(m, e1)
                return true
            end
        end
    end

    return false

end

function p2_bisect_edge!(m::Mesh, edge::Base.RefValue{Edge})

    # Generate new node
    new_node = Ref(Node(new_coords(edge), new_coords(edge)))

    push!(m.nodes, new_node)
    new_node_id = Base.length(m.nodes)

    # Generate the first new edge
    edge1 = Ref(Edge([edge.x.nodes[1], new_node], [edge.x.nodes_id[1], new_node_id], false, edge.x.NA, nothing))

    # Generate the second new edge
    edge2 = Ref(Edge([new_node, edge.x.nodes[2]], [new_node_id, edge.x.nodes_id[2]], false, edge.x.NA, nothing))

    # Add the sons of the initial edge
    edge.x.sons = [edge1, edge2]

    return nothing

end

function p3_bisect_triangle!(m::Mesh, triangle::Base.RefValue{Triangle})

    if isbroken(triangle.x)
        return false
    end

    edge1, edge2, edge3 = get_sorted_edges(triangle.x)

    if isbroken(edge1)
        triangle.x.MR = false
        triangle.x.BR = true

        edge4, edge5 = edge1.x.sons

        if (isempty(common_node(edge3, edge4)))
            edge4, edge5 = edge5, edge4
        end

        v1 = common_node(edge4, edge5)[]
        v2 = common_node(edge2, edge3)[]
        v1_id = common_node_id(edge4, edge5)[]
        v2_id = common_node_id(edge2, edge3)[]
        new_edge = Ref(Edge([v1, v2], [v1_id, v2_id], false, 2, nothing))

        triangle1 = Triangle([edge3, edge4, new_edge], false, false, nothing, nothing)
        triangle2 = Triangle([edge5, edge2, new_edge], false, false, nothing, nothing)

        add_new_triangles!(m, triangle.x, triangle1, triangle2)

        return true
    else
        return false
    end

end

function collect_all_triangles(m::Mesh)
    node = m.root_triangle
    triangles = Vector{Base.RefValue{Triangle}}()
    while !isnothing(node)
        push!(triangles, node)
        node = node.x.next
    end
    return triangles
end

function refine!(m::Mesh)

    run = true

    while run
        run = false
        l_triangles = collect_all_triangles(m)

        Threads.@threads for t in l_triangles
            run |= p1_mark_edges!(m, t)
            run |= p3_bisect_triangle!(m, t)
        end

    end

    return nothing

end

collect_all_nodes(m::Mesh) = m.nodes
get_conec(t::Triangle) = [n for n in 
    [
        common_node_id(t.edges[3], t.edges[1])[],
        common_node_id(t.edges[1], t.edges[2])[],
        common_node_id(t.edges[2], t.edges[3])[]
    ]
]

function write_vtk(m::Mesh, filename)
    nodes_vec = collect_all_nodes(m)
    coords = mapreduce(n -> n.x.xyz, hcat, nodes_vec)
    conec = [MeshCell(VTKCellTypes.VTK_TRIANGLE, get_conec(t.x)) for t in collect_all_triangles(m)]
    vtk_grid(filename, coords, conec) do vtk
        vtk["uvw"] = mapreduce(n -> n.x.uvw, hcat, nodes_vec)
    end
end

end
