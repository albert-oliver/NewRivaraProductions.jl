module NewRivaraProductions

using LinearAlgebra

# Write your package code here.

mutable struct Node
    uvw::Vector{Float64}
    xyz::Vector{Float64}
end

mutable struct Edge
    v1::Base.RefValue{Node}
    v2::Base.RefValue{Node}
    MR::Bool
    BR::Bool
    NA::Int
    sons::Vector{Int}
end

mutable struct Triangle
    edges::Vector{Int}
    MR::Bool
    BR::Bool
end

mutable struct Mesh
    nodes::Vector{Node}
    all_edges::Vector{Edge}
    all_triangles::Vector{Triangle}
    edges::Vector{Int}
    triangles::Vector{Int}
    function Mesh(n::Vector{Node}, a_t::Vector{Triangle})
        m = new()
        m.nodes = n
        m.all_triangles = a_t
        m.triangles = 1:Base.length(a_t)
        m.all_edges = []
        m.edges = []
        return m
    end
end

function simpleMesh()
    m = Mesh([Node([0, 0, 0], [0, 0, 0]), Node([1, 0, 0], [1, 0, 0]),
        Node([1, 1, 0], [1, 1, 0]), Node([0, 1, 0], [0, 1, 0])],
    [Triangle([1, 2, 5], false, false), Triangle([3, 4, 5], false, false)]
    )

    n1 = Ref(m.nodes[1])
    n2 = Ref(m.nodes[2])
    n3 = Ref(m.nodes[3])
    n4 = Ref(m.nodes[4])
    push!(m.all_edges, Edge(n1, n2, false, false, 1, []))
    push!(m.all_edges, Edge(n2, n3, false, false, 1, []))
    push!(m.all_edges, Edge(n3, n4, false, false, 1, []))
    push!(m.all_edges, Edge(n4, n1, false, false, 1, []))
    push!(m.all_edges, Edge(n1, n3, false, false, 2, []))

    m.edges = 1:Base.length(m.all_edges)
    return m
end

length(e::Edge) = norm(e.v1.x.xyz - e.v2.x.xyz)
new_coords(e::Edge) = (e.v1.x.xyz + e.v2.x.xyz)/2

common_node(e1::Edge, e2::Edge) = intersect([e1.v1, e1.v2], [e2.v1, e2.v2])

istrianglebroken(m::Mesh, t::Int) = m.all_triangles[t].BR
istrianglebroken(m::Mesh) = Base.Fix1(istrianglebroken, m)

isedgebroken(m::Mesh, e::Int) = m.all_edges[e].BR
isedgebroken(m::Mesh) = Base.Fix1(isedgebroken, m)

function Base.isless(e1::Edge, e2::Edge)

    # First we just compare lengths
    l1 = length(e1)
    l2 = length(e2)
    if l1 ≉ l2
        return l1 < l2
    end

    # Then, if any edge is marked to be refined (or already broken)
    mr1 = (e1.MR || e1.BR)
    mr2 = (e2.MR || e2.BR)
    if mr1 ≠ mr2
        # This function returns the edge that is NOT going to be refined,
        # we want to refine the edge that is marked for refinement or already broken,
        # therefore if it's marked or already broken it should return false
        return !(mr1)
    end

    # Next, the number of adjacent elements
    na1 = e1.NA
    na2 = e2.NA
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

function get_edge_index(m::Mesh, e::Edge)
    return findfirst(isequal(e), m.all_edges)
end

function get_edge(m::Mesh, e::Int)
    return m.all_edges[e]
end

function get_edge(m::Mesh, e::Vector{Int})
    return m.all_edges[e]
end

function get_triangle(m::Mesh, t::Int)
    return m.all_triangles[t]
end

function get_triangle_edges(m::Mesh, t::Int)
    tri = m.all_triangles[t]
    return m.all_edges[tri.edges]
end

function add_new_node!(m::Mesh, new_node::Node)
    push!(m.nodes, new_node)
    return Ref(m.nodes[findlast(isequal(new_node), m.nodes)])
end

function add_new_edge!(m::Mesh, e::Edge)
    push!(m.all_edges, e)
    e1 = findlast(isequal(e), m.all_edges)
    push!(m.edges, e1)
    return e1
end

function add_new_triangle!(m::Mesh, t::Triangle)
    push!(m.all_triangles, t)
    t1 = findlast(isequal(t), m.all_triangles)
    push!(m.triangles, t1)
    return t1
end

function p1_mark_edges!(m::Mesh, t::Int)
    triangle = get_triangle(m, t)

    if (triangle.BR)
        return false
    end

    e1, e2, e3 = sort(get_triangle_edges(m, t), rev=true)

    if (!e1.MR) && (triangle.MR || (e2.MR || e2.BR) || (e3.MR || e3.BR))
        e1.MR = true
        return true
    else
        return false
    end

end

function p2_bisect_edge!(m::Mesh, e::Int)
    edge = get_edge(m, e)

    if edge.MR
        # It's not marked for refinement anymore
        edge.MR = false
        # because it's broken
        edge.BR = true

        # Generate new node
        new_node = Node(new_coords(edge), new_coords(edge))
        v3 = add_new_node!(m, new_node)

        # Generate the first new edge
        edge1 = Edge(edge.v1, v3, false, false, edge.NA, [])
        e1 = add_new_edge!(m, edge1)

        # Generate the second new edge
        edge2 = Edge(v3, edge.v2, false, false, edge.NA, [])
        e2 = add_new_edge!(m, edge2)

        # Add the sons of the initial edge
        edge.sons = [e1, e2]

        # Production has been performed
        return true
    else
        # This production has not been performed
        return false
    end
end

function p3_bisect_triangle!(m::Mesh, t::Int)
    triangle = get_triangle(m, t)

    if (triangle.BR)
        return false
    end

    edge1, edge2, edge3 = sort(get_triangle_edges(m, t), rev=true)


    if edge1.BR
        triangle.MR = false
        triangle.BR = true

        edge4, edge5 = get_edge(m, edge1.sons)

        if (isempty(common_node(edge3, edge4)))
            edge4, edge5 = edge5, edge4
        end

        e2 = get_edge_index(m, edge2)
        e3 = get_edge_index(m, edge3)
        e4 = get_edge_index(m, edge4)
        e5 = get_edge_index(m, edge5)

        v1 = common_node(edge4, edge5)[]
        v2 = common_node(edge2, edge3)[]
        new_edge = Edge(v1, v2, false, false, 2, [])
        e6 = add_new_edge!(m, new_edge)

        triangle1 = Triangle([e3, e4, e6], false, false)
        add_new_triangle!(m, triangle1)

        triangle2 = Triangle([e5, e2, e6], false, false)
        add_new_triangle!(m, triangle2)

        return true

    else
        return false
    end

end

p4_remove_broken_edges!(m::Mesh) = filter!(!isedgebroken(m), m.edges)

p5_remove_broken_triangles!(m::Mesh) = filter!(!istrianglebroken(m), m.triangles)

function refine!(m::Mesh)

    run = true

    while run
        run = false

        for t in m.triangles
            run |= p1_mark_edges!(m, t)
        end

        for e in m.edges
            run |= p2_bisect_edge!(m, e)
        end

        for t in m.triangles
            run |= p3_bisect_triangle!(m, t)
        end

    end

    p4_remove_broken_edges!(m)
    p5_remove_broken_triangles!(m)

end

end
