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
    MR::Bool
    NA::Int
    sons::Union{SVector{2,Base.RefValue{Edge}}, Nothing}
end

mutable struct Triangle
    edges::SVector{3,Base.RefValue{Edge}}
    MR::Bool
    BR::Bool
end

mutable struct Mesh
    triangles::Vector{Triangle}
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
                edge = Ref(Edge([nodes[n1], nodes[n2]], false, 1, nothing))
                push!(get!(edges_per_node, n1, []), edge)
                push!(get!(edges_per_node, n2, []), edge)
            else
                edge = edge[1]
                edge.x.NA = edge.x.NA + 1
            end
            edges[j] = edge
        end
        triangles[i] = Triangle(edges, false, false)
    end
    return Mesh(triangles)
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

function get_triangle_edges(t::Triangle)
    return t.edges
end

lk = ReentrantLock()

function add_new_triangle!(m::Mesh, t::Triangle)
    lock(lk) do
        push!(m.triangles, t)
    end
end

function p1_mark_edges!(triangle::Triangle)

    if isbroken(triangle)
        return false
    end

    e1, e2, e3 = sort(get_triangle_edges(triangle), rev=true)

    if (!isbroken(e1)) && (triangle.MR || e1.x.MR || (isbroken(e2) || e2.x.MR) || (isbroken(e3) || e3.x.MR))
        return p2_bisect_edge!(e1)
    else
        return false
    end

end

function p2_bisect_edge!(edge::Base.RefValue{Edge})

    # Generate new node
    new_node = Ref(Node(new_coords(edge), new_coords(edge)))

    # Generate the first new edge
    edge1 = Ref(Edge([edge.x.nodes[1], new_node], false, edge.x.NA, nothing))

    # Generate the second new edge
    edge2 = Ref(Edge([new_node, edge.x.nodes[2]], false, edge.x.NA, nothing))

    # Add the sons of the initial edge
    lock(lk) do
        if !isbroken(edge)
            edge.x.sons = [edge1, edge2]
            return true
        end
    end

    return false

end

function p3_bisect_triangle!(m::Mesh, triangle::Triangle)

    if isbroken(triangle)
        return false
    end

    edge1, edge2, edge3 = sort(get_triangle_edges(triangle), rev=true)


    if isbroken(edge1)
        triangle.MR = false
        triangle.BR = true

        edge4, edge5 = edge1.x.sons

        if (isempty(common_node(edge3, edge4)))
            edge4, edge5 = edge5, edge4
        end

        v1 = common_node(edge4, edge5)[]
        v2 = common_node(edge2, edge3)[]
        new_edge = Ref(Edge([v1, v2], false, 2, nothing))

        triangle1 = Triangle([edge3, edge4, new_edge], false, false)
        add_new_triangle!(m, triangle1)

        triangle2 = Triangle([edge5, edge2, new_edge], false, false)
        add_new_triangle!(m, triangle2)

        return true

    else
        return false
    end

end

p4_remove_broken_triangles!(m::Mesh) = filter!(!isbroken, m.triangles)

function collect_all_edges(m)
    return collect(Set(reduce(vcat, [t.edges for t in m.triangles])))
end

function refine!(m::Mesh)

    run = true

    while run
        run = false
        l_triangles = copy(m.triangles)

        Threads.@threads for t in l_triangles
            run |= p1_mark_edges!(t)
            run |= p3_bisect_triangle!(m, t)
        end

    end

    p4_remove_broken_triangles!(m)
    return nothing

end

collect_all_nodes(m::Mesh) = collect_all_nodes(collect_all_edges(m))
collect_all_nodes(edge_vector::AbstractVector{Base.RefValue{Edge}}) = collect(Set(reduce(vcat, [e.x.nodes for e in edge_vector])))
node_idx_Dict(node_vector::Vector{Base.RefValue{Node}}) = Dict([(node, i) for (i, node) in enumerate(node_vector)])
get_triangle_conec(t::Triangle, nodes_dict) = [nodes_dict[n] for n in collect_all_nodes(t.edges)]

function write_vtk(m::Mesh, filename)
    nodes_vec = collect_all_nodes(m)
    nodes_dict = node_idx_Dict(nodes_vec)
    coords = mapreduce(n -> n.x.xyz, hcat, nodes_vec)
    conec = [MeshCell(VTKCellTypes.VTK_TRIANGLE, get_triangle_conec(t, nodes_dict)) for t in m.triangles]
    vtk_grid(filename, coords, conec) do vtk
        vtk["uvw"] = mapreduce(n -> n.x.uvw, hcat, nodes_vec)
    end
end

end
