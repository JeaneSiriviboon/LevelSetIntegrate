
using Plots;
using Meshes;
using StatsBase;
import CairoMakie as Mke;
using SpecialFunctions;
using BenchmarkTools;
using MultiQuad;

"""
Data Structure 

Denote grid with tree

Mesh: List of mesh and unit mesh
UnitMesh: triangular mesh
"""
abstract type Mesh end
abstract type Vertex end
abstract type Surface end
# abstract type UnitMesh2D          <: Mesh2D end

struct Triangle 
    v1::Int
    v2::Int 
    v3::Int
end

struct Vertex2D <: Vertex
    x::Float64
    y::Float64
end

struct Vertex2DComplex <: Vertex
    x::ComplexF64
    y::ComplexF64
end


Base.:+(v1::Vertex2D, v2::Vertex2D)::Vertex2D = Vertex2D(v1.x + v2.x, v1.y + v2.y)
Base.:+(c::Float64, v2::Vertex2D)::Vertex2D = Vertex2D(c + v2.x, c + v2.y)
Base.:-(v1::Vertex2D, v2::Vertex2D)::Vertex2D = Vertex2D(v1.x - v2.x, v1.y - v2.y)
Base.:/(v::Vertex2D, c::Float64)::Vertex2D = Vertex2D(v.x /c, v.y /c)
Base.:*(c::Float64, v::Vertex2D)::Vertex2D = Vertex2D(v.x *c, v.y*c)
Base.:*(v::Vertex2D, c::Float64)::Vertex2D = Vertex2D(v.x *c, v.y*c)

Base.:+(v1::Vertex2DComplex, v2::ComplexF64)::Vertex2DComplex = Vertex2DComplex(v1.x + v2.x, v1.y + v2.y)
Base.:+(c::ComplexF64, v2::Vertex2DComplex)::Vertex2DComplex = Vertex2DComplex(c + v2.x, c + v2.y)
Base.:-(v1::Vertex2DComplex, v2::Vertex2DComplex)::Vertex2DComplex = Vertex2DComplex(v1.x - v2.x, v1.y - v2.y)
Base.:/(v::Vertex2DComplex, c::ComplexF64)::Vertex2DComplex = Vertex2DComplex(v.x /c, v.y /c)
Base.:*(c::ComplexF64, v::Vertex2DComplex)::Vertex2DComplex = Vertex2DComplex(v.x *c, v.y*c)
Base.:*(v::Vertex2DComplex, c::ComplexF64)::Vertex2DComplex = Vertex2DComplex(v.x *c, v.y*c)

struct Edge
    p1::Vertex
    p2::Vertex
end

struct Mesh2D <: Mesh
    _vertices::Vector{Vertex2D}
    _triangles::Vector{Triangle}
    
    function Mesh2D()
        new([], [])
    end 

    function Mesh2D(_vertices, _triangles)
        new(_vertices, _triangles)
    end 
end

struct Surface2D<:Surface
    _mesh::Mesh2D 
    _values::Vector{Float64}
    _Q_fn::Any

    function Surface2D(mesh::Mesh2D, f::Any)
        values = []
        for pos in mesh._vertices
            Base.push!(values, f(pos))
        end
        new(mesh, values, f)
    end
end



function Base.copy(mesh::Mesh2D)
    return Mesh2D(copy(mesh._vertices), copy(mesh._triangles))
end

function get_area(mesh::Mesh2D, i::Int)::Float64
    v1 = mesh._triangles[i].v1 
    v2 = mesh._triangles[i].v2 
    v3 = mesh._triangles[i].v3
    
    u1 = mesh._vertices[v2] - mesh._vertices[v1]
    u2 = mesh._vertices[v3] - mesh._vertices[v1]

    area = (u1.x * u2.y - u1.y * u2.x)/2
    return area
end

function get_square_grid(X, Y, n, rand_pos = 0.0)::Mesh2D
    dx = X/n
    dy = Y/n
    mesh = Mesh2D()

    #(i, j) vertices locate at (i - 1) * (n+1) + j 
    for i = 1:n+1 
        for j = 1:n+1
            x = (i-1)*dx
            y = (j-1)*dy
            
            if i > 1 && i < n+1
                x = clamp(x * (1 - rand_pos/2 + rand_pos*rand()), 0, X)
            end
            
            if j > 1 && j < n+1
                y = clamp(y * (1 - rand_pos/2 + rand_pos*rand()), 0, Y)
            end

            Base.push!(mesh._vertices, 
                        Vertex2D(x, y))

            if i < n+1  && j < n+1
                v1 = (i - 1) * (n+1) + j
                v2 = (i) * (n+1) + j
                v3 = (i - 1) * (n+1) + (j+1)
                v4 = (i) * (n+1) + (j+1)
                v5 = (n+1)*(n+1) + (i - 1) * n + j

                Base.push!(mesh._triangles,Triangle(v1, v2, v5))  
                Base.push!(mesh._triangles,Triangle(v1, v5, v3))
                Base.push!(mesh._triangles,Triangle(v2, v4, v5))
                Base.push!(mesh._triangles,Triangle(v4, v3, v5))
            end
        end
    end

    for i=1:n
        for j=1:n
            x = (i-0.5)*dx* (1 - rand_pos/2 + rand_pos*rand())
            y = (j-0.5)*dy* (1 - rand_pos/2 + rand_pos*rand())
            Base.push!(mesh._vertices, 
                        Vertex2D(x, y))
        end
    end

    return mesh
end

function plot(mesh::Mesh2D)
    points = [(p.x, p.y) for p in mesh._vertices]
    connec = [connect((t.v1 , t.v2, t.v3)) for t in mesh._triangles]
    simple_mesh   = SimpleMesh(points, connec)
    Meshes.viz(simple_mesh, showsegments = true)
end


function plot(surf::Surface2D)
    points = [(p.x, p.y) for p in surf._mesh._vertices]
    connec = [connect((t.v1 , t.v2, t.v3)) for t in surf._mesh._triangles]
    colors = [surf._values[i] for i = 1:length(surf._mesh._vertices)]
    simple_mesh   = SimpleMesh(points, connec)
    # println(simple_mesh)
    vmax = StatsBase.percentile(abs.(colors), 95)
    Meshes.viz(simple_mesh, 
                color = clamp.(colors, -vmax, vmax), 
                showsegments = true, 
                segmentcolor = RGBA(0,0,0,0.1),
                colormap = :balance)
end

function subsample!(mesh::Mesh2D, i::Int)
    t = mesh._triangles[i]
    v1 = mesh._vertices[t.v1]
    v2 = mesh._vertices[t.v2]
    v3 = mesh._vertices[t.v3]
    
    u1 = Vertex2D((v2.x + v3.x) /2, 
                  (v2.y + v3.y) /2)
    u2 = Vertex2D((v3.x + v1.x) /2, 
                  (v3.y + v1.y) /2)
    u3 = Vertex2D((v1.x + v2.x) /2, 
                  (v1.y + v2.y) /2)

    idx1 = length(mesh._vertices) + 1
    idx2 = length(mesh._vertices) + 2
    idx3 = length(mesh._vertices) + 3

    Base.push!(mesh._vertices, u1)
    Base.push!(mesh._vertices, u2)
    Base.push!(mesh._vertices, u3)

    t0 = Triangle(idx1, idx2, idx3)
    t1 = Triangle(t.v1, idx3, idx2)
    t2 = Triangle(t.v2, idx1, idx3)
    t3 = Triangle(t.v3, idx2, idx1)
    # print(t0)
    
    Base.push!(mesh._triangles, t1)
    Base.push!(mesh._triangles, t2)
    Base.push!(mesh._triangles, t3)
    mesh._triangles[i] = t0
end

function subsample!(surf::Surface2D, i::Int)
    L = length(surf._mesh._vertices)

    subsample!(surf._mesh, i)
    
    u1, u2, u3 = surf._mesh._vertices[L+1], surf._mesh._vertices[L+2], surf._mesh._vertices[L + 3]
    
    Base.push!(surf._values, surf._Q_fn(u1))
    Base.push!(surf._values, surf._Q_fn(u2))
    Base.push!(surf._values, surf._Q_fn(u3))
end


function subsample_check(surf::Surface2D, abs_err::Float64, i::Int)::Bool
    m = surf._mesh._triangles[i]

    v1, v2, v3 = mesh._vertices[m.v1], mesh._vertices[m.v2], mesh._vertices[m.v3]

    integrand1, integrand2, integrand3 = surf._values[m.v1], surf._values[m.v2], surf._values[m.v3]

    x_center = (v1.x + v2.x + v3.x)/3
    y_center = (v1.x + v2.x + v3.x)/3

    integrand4 = surf._Q_fn(x_center, y_center)      
    
    l1_sq = (v2.x -v1.x)^2 + (v2.y -v1.y)^2
    l2_sq = (v3.x -v1.x)^2 + (v3.y -v1.y)^2

    area = (l2_sq * l1_sq)^0.5

    err = area * abs(integrand4 - (integrand1 + integrand2 + integrand3)/3)

    success = err < abs_err

    return success
end

function subsample!(surf::Surface2D, abs_err::Float64, depth::Int)
    idxs::Vector{Int} = 1:length(surf._mesh._triangles)
    for n = 1:depth
        idx_new::Vector{Int} = []
        if length(idxs) > 0
            for i in idxs
                success = subsample_check(surf, abs_err, i)
                if !success
                    subsample!(surf, i)
                    L = length(surf._mesh._triangles)
                    push!(idx_new, i)
                    push!(idx_new, L-1)
                    push!(idx_new, L-2)
                    push!(idx_new, L-3)
                end
            end
        end
        idxs = idx_new
    end
end


function interp(u::Vertex2D, v::Vertex2D, p::Float64)
    return Vertex2D(p * u.x + (1-p) * v.x, p * u.y + (1-p) * v.y)
end

function interp(u::Vertex2DComplex, v::Vertex2DComplex, p::Float64)
    return Vertex2DComplex(p * u.x + (1-p) * v.x, p * u.y + (1-p) * v.y)
end


norm(v::Vertex2D) = sqrt(v.x^2 + v.y^2)
norm(v::Vertex2DComplex) = sqrt(abs(v.x)^2 + abs(v.y)^2)


function num_singular(vals::Vector{Float64})
    """
    check number of zeros within edges and verts
    """
    num_zeros_vert = Int(vals[1] == 0.) + Int(vals[2] == 0.) + Int(vals[3] == 0.)

    num_zeros_edge = Int(vals[1] * vals[2] < 0.) + 
                     Int(vals[2] * vals[3] < 0.) + 
                     Int(vals[1] * vals[3] < 0.)

    return num_zeros_edge + num_zeros_vert 
end

function num_singular(surf::Surface2D, c::Float64, i::Int, tol::Float64)
    t = surf._mesh._triangles[i]
    fs  = [surf._values[t.v1], surf._values[t.v2], surf._values[t.v3]]
    vals = fs .- c 
    return num_singular(vals) 
end


function is_subsample(surf::Surface2D, c::Float64, i::Int, tol::Float64)
    """
    check if the error ε =  | Q((vi + vj)/2) - (Q(vi) + Q(vj))/2 | exceeds the threshold
    """
    t = surf._mesh._triangles[i]
    fs  = [surf._values[t.v1], surf._values[t.v2], surf._values[t.v3]]
    ps  = [surf._mesh._vertices[t.v1], surf._mesh._vertices[t.v2], surf._mesh._vertices[t.v3]]
    vals = fs .- c 
    #value average for each edges
    vals_avg = [(vals[2] + vals[3])/2., (vals[3] + vals[1])/2., (vals[1] + vals[2])/2.]
    
    pc = [(ps[2] + ps[3])/2., (ps[3] + ps[1])/2., (ps[1] + ps[2])/2.]

    vals_c = surf._Q_fn.(pc) .- c
    # dl = [norm(ps[2] - ps[3]), norm(ps[3] - ps[1]), norm(ps[1] - ps[2])]
    
    err = abs.(vals_avg .- vals_c)

    return any(err .> tol)
end

function level_set_subsample!(surf::Surface, 
                                c::Float64, 
                                depth::Int, 
                                tol::Float64)::Vector{Int}
    """
    Subsampling the mesh until we get the level set
    """
    idxs = 1:length(surf._mesh._triangles) 
    idx_keep = []
    idx_new = []
    for n =1:depth - 1
        idx_new = []
        for i in idxs

            sample = is_subsample(surf, c, i, tol)
            sing =   num_singular(surf, c, i, tol) > 0

            if sample && sing
                push!(idx_new, i)
                subsample!(surf, i)
                L = length(surf._mesh._triangles)
                push!(idx_new, L-3)
                push!(idx_new, L-2)
                push!(idx_new, L-1)
            elseif sing
                push!(idx_keep, i)
            end

        end
        idxs = idx_new
    end
    idxs = Base.vcat(idx_keep, idx_new)
    return idxs
end



function get_gradient(surf::Surface, i::Int)
    """
    Finite element method for calculating |∇Q|
    """
    t = surf._mesh._triangles[i]
    fs  = [surf._values[t.v1], surf._values[t.v2], surf._values[t.v3]]
    ps  = [surf._mesh._vertices[t.v1], surf._mesh._vertices[t.v2], surf._mesh._vertices[t.v3]]
    df  = [fs[2] - fs[1], fs[3] - fs[1]]
    dp  = [ps[2] - ps[1], ps[3] - ps[1]]

    dp1dp2 = dp[1].x * dp[2].y -  dp[2].x * dp[1].y

    gx = dp[2].y * df[1] - dp[1].y * df[2]
    gy = - dp[2].x * df[1] + dp[1].x * df[2]

    return sqrt(gx^2 + gy^2)/dp1dp2
end

function level_set(surf::Surface, c::Float64, i::Int)
    """
    Calculated list of edge that forms a level set for each triangle
    """
    t = surf._mesh._triangles[i]
    fs  = [surf._values[t.v1], surf._values[t.v2], surf._values[t.v3]]
    ps  = [surf._mesh._vertices[t.v1], surf._mesh._vertices[t.v2], surf._mesh._vertices[t.v3]]
    vs = fs .- c

    zeros = []

    for i = 1:3
        if vs[i] == 0. 
            push!(zeros, ps[i])
        end

        j = i%3 + 1
        if vs[i] * vs[j] < 0. 
            p = vs[j]/(vs[j] - vs[i])
            zero_root = interp(ps[i], ps[j], p)
            push!(zeros, zero_root)
        end
    end

    edges::Vector{Edge} = []

    L = length(zeros)
    
    for i = 1: L
        if i!=2 || L!=2 
            push!(edges, Edge(zeros[i], zeros[i %L + 1]))
        end
    end
    return edges
end


function level_set(surf::Surface, c::Float64, idxs::Vector{Int}, tol::Float64)::Vector{Edge}
    """
    Calculated list of edge that forms a level set for the whole mesh
    """
    level_set_edges = []
        
    for i in idxs
        sing = num_singular(surf, c, i, tol) >= 2
        if sing  
            edges = level_set(surf, c, i)
            for edge in edges
                push!(level_set_edges, edge)
            end
        end
    end

    return level_set_edges
end

function level_set!(surf::Surface, c::Float64, depth::Int, tol::Float64)::Vector{Edge}
    """
    Calculated list of edge that forms a level set for the whole mesh
    """
    idxs = level_set_subsample!(surf, c, depth, tol)

    level_set_edges = level_set(surf, c, idxs, tol)

    return level_set_edges
end

conj(v::Vertex2DComplex) = Vertex2DComplex(Base.conj(v.x), Base.conj(v.y))

#TODO: Write early stopping
function find_zeros(p::Vertex2D, f::Any, grad::Any, η::Any ,depth::Int)::Vertex2DComplex
    # p1 s.t. find f(p1) + iη(p1) = 0.0 + 0.0im
    p1 = Vertex2DComplex(ComplexF64(p.x), ComplexF64(p.y))
    for i = 1:depth
        
        grad_vec = grad(p1) 
        p1 -= (f(p1) + 1.0im *η(p1))/norm(grad_vec)^2 * conj(grad_vec) 
    end

    return p1
end

function plot_level_set!(l_set::Vector{Edge})
    
    for edge in l_set 
        segment = Meshes.Segment((edge.p1.x, edge.p1.y), 
                                (edge.p2.x, edge.p2.y))
        viz!(segment, color = RGBA(50,0,50,0.5))
    end
end

function integrate1D_helper(f::Any, a::Float64, b::Float64)
    # println(a, " ", b)
    fl = f(a)
    # println("f(a)")
    fr = f(b)
    fc = f((a + b)/2.)
    L = b-a
    err = abs(fc - (fl + fr)/2.)

    res = (fl + 4 * fc + fr)/6 * L

    return res, err
end

function integrate1D(a::Float64, b::Float64, 
                     f_fn::Any, depth::Int, tol::Float64)
   

    xl = [a]
    xr = [b]
    acc = 0.

    for n=1:depth 
        xl_ = []
        xr_ = []
        Len = length(xl)
        
        for i = 1:Len

            # println("xlxr:", xl[i], " ", xr[i])
            res, err = integrate1D_helper(f_fn, xl[i], xr[i])
            # println("res:", res)
            
            if (err < tol || n == depth) && !isnan(res)
                acc += res
            else
                push!(xl_, xl[i])
                push!(xr_, (xl[i] + xr[i])/2.)
                
                push!(xl_, (xl[i] + xr[i])/2.)
                push!(xr_, xr[i])
            end
        end

        xl = xl_
        xr = xr_
    end

    return acc
end


function integrate1D(e::Edge, f_fn::Any, depth::Int, tol::Float64)
    vl, vr = e.p1, e.p2
    L = norm(vr - vl)

    x2p(p::Float64) = interp(vl, vr, 1 - p)
    integrand(p::Float64) = f_fn(x2p(p))
    return L * integrate1D(0., 1., integrand, depth, tol)
end

# function integrate1D(e::Edge, grad::Vector{Float64}, f::Any, depth::Int, tol::Float64)
#     vl, vr = e.p1, e.p2
    
#     L = norm(vr - vl)

#     x2p(p) = interp(vl, vr, 1 - p)

#     f_grad(x) = f(x2p(x))/(grad[1] + x/L * grad[2])

#     return integrate1D(e, f_grad, depth, tol)
# end

function integrate1D(e::Edge, grad::Any, f::Any, depth::Int, tol::Float64)
    # println("This should be called")
    vl, vr = e.p1, e.p2

    norm_analytic(g::Vertex)::Union{ComplexF64, Float64} = sqrt((g.x)^2 + (g.y)^2) 

    f_grad(x) = f(x)/norm_analytic(grad(x))

    return integrate1D(e, f_grad, depth, tol)
end


function integrate1D(l_set::Vector{Edge}, grad::Vector{Vector{Float64}}, p::Any, depth::Int, tol::Float64)
    acc = 0

    for i = 1:length(l_set)
        acc += integrate1D(l_set[i], grad[i], p, depth,  tol)
    end
    return acc
end


function integrate_delta(p_fn::Any, q_fn::Any, domain::Mesh2D; 
                        depth_LS  = 6,
                        depth_int = 6,
                        tol_LS = 1e-9, 
                        tol_int = 1e-9)
    """
    Calculate ∫dx P(x)δ(Q(x))
    """

    surf = Surface2D(copy(domain), q_fn)
    level_set_subsample!(surf, 0., depth_LS, tol_LS)
    L = length(surf._mesh._triangles)
    res = 0.

    # looping for each triangle
    for i = 1:L 
        # check for isoline
        if num_singular(surf, 0., i, tol_LS) >= 2
            grad = get_gradient(surf, i)   # get gradient for the triangle
            edges = level_set(surf, 0., i) # get the edge
            for edge in edges
                res += integrate1D(edge, p_fn, depth_int, tol_int)/grad
            end
        end     
    end
    return res 
end

function integrate_im(p_fn::Any, q_fn::Any, Domain::Mesh2D; 
    depth_LS  = 7,
    depth_int = 7,
    tol_LS = 1e-10, 
    tol_int = 1e-10)
    return -π * integrate_delta(p_fn, q_fn, Domain; depth_LS, depth_int, tol_LS, tol_int)
end

function integrate_non_sing_helper(f1::Float64, 
                                   f2::Float64, 
                                   f3::Float64, 
                                   fc::Float64, 
                                   area::Float64)
        
        err = abs(fc - (f1 + f2 + f3)/3.)
        res = area * ((f1 + f2 + f3) * 2. / 9. + fc/3.)
        return res, err
end

function integrate_non_sing(surf_Q::Surface2D, 
                            p_fn::Any, 
                            idxs::Vector{Int}, 
                            depth::Int, 
                            tol::Float64)
    mesh_Q = surf_Q._mesh
    
    p_vals = p_fn.(mesh_Q._vertices)

    acc = 0

    for n = 1:depth
        idx_new = []
        for j in idxs

            t = mesh_Q._triangles[j]
            
            x1, x2, x3 = mesh_Q._vertices[t.v1], 
                         mesh_Q._vertices[t.v2], 
                         mesh_Q._vertices[t.v3]

            p1, p2, p3 = p_vals[t.v1], p_vals[t.v2], p_vals[t.v3]
            
            q1, q2, q3 = surf_Q._values[t.v1],
                         surf_Q._values[t.v2],
                         surf_Q._values[t.v3]

            f1, f2, f3 = p1/q1, p2/q2, p3/q3
            xc = (x1 + x2 + x3)/3.
            fc = p_fn(xc)/surf_Q._Q_fn(xc)

            area = get_area(mesh_Q, j)
            
            res, err = integrate_non_sing_helper(f1, f2, f3, fc, area)
            if n == depth || err < tol
                acc += res
            else
                Lt = length(surf_Q._mesh._triangles)
                Lv = length(surf_Q._mesh._vertices)
                subsample!(surf_Q, j)
                push!(idx_new, j)
                push!(idx_new, Lt - 3)
                push!(idx_new, Lt - 2)
                push!(idx_new, Lt - 1)
                push!(p_vals, p_fn(mesh_Q._vertices[Lv - 3]))
                push!(p_vals, p_fn(mesh_Q._vertices[Lv - 2]))
                push!(p_vals, p_fn(mesh_Q._vertices[Lv - 1]))
            end
        end
        idxs = idx_new
    end

    return acc
end

function integrate_re(p_fn::Any, q_fn::Any, domain::Mesh2D; 
    depth_LS  = 4,
    depth_int = 4,
    tol_LS = 1e-6, 
    tol_int = 1e-6,
    lnc_min = -6.,
    dx = 1,
    tol_inf = 1e-7)
    
    surf = Surface2D(domain, q_fn)
    lnc_max = log(maximum(abs.(surf._values))) + dx


    function integrand(lnc)
        c = exp(lnc)
        int_pos = integrate_delta(p_fn, x->(q_fn(x) - c), 
                                  domain; 
                                  depth_LS, 
                                  depth_int, 
                                  tol_LS, 
                                  tol_int)
        int_neg = integrate_delta(p_fn, x->(q_fn(x) + c), 
                                  domain; 
                                  depth_LS, 
                                  depth_int, 
                                  tol_LS, 
                                  tol_int)
        return int_pos - int_neg
    end

    acc = 0.
    threshold = -1.
    for lnc = lnc_max:-dx:lnc_min
        
        res = integrate1D(lnc, lnc + dx, integrand, depth_int, tol_int)    
        acc += res
        if abs(res)/max(abs(acc), 1e-10) < tol_inf && lnc < threshold
            break
        end
    end

    return acc
end

function calculate_dos(ω::Float64, ε_k::Any, BZ::Mesh2D; 
    depth_LS  = 7,
    depth_int = 7,
    tol_LS = 1e-10, 
    tol_int = 1e-10)::Float64 
    return  integrate_delta(k -> 1., k -> (ω - ε_k(k)), BZ; depth_LS, depth_int, tol_LS, tol_int)
end



function mesh_test1()

    function p_fn(v)
        return 1.0
    end

    function q_fn(v) 
        return  cos.(2* π *  v.x) + cos.(2 * π * v.y)
    end

    function f_fn(v)
        return q_fn(v)
    end

    t0  = time()
    mesh  = get_square_grid(1., 1., 31)
    surf = Surface2D(mesh, q_fn)

    p = plot(surf)

    cmax = max(maximum(surf._values), 0)
    cmin = min(minimum(surf._values), 0)

    if cmax > 0
        for lnc = -6.:1.: log(cmax)+1.
            mesh  = get_square_grid(1., 1., 57)
            surf = Surface2D(mesh, q_fn)
            level_set_edges = level_set!(surf, exp(lnc), 4, 1e-4)
            plot_level_set!(level_set_edges)
            println("time: ", time() - t0," lnc: ", lnc, " level set: ", length(level_set_edges))
        end
    end


    if cmin < 0 
        for lnc = -6.:1.: log(-cmin)+1.
            mesh  = get_square_grid(1., 1., 57)
            surf = Surface2D(mesh, q_fn)
            level_set_edges = level_set!(surf, -exp(lnc), 4, 1e-4)
            plot_level_set!(level_set_edges)
            println("time: ", time() - t0," lnc: ", lnc, " level set: ", length(level_set_edges))
        end
    end
    Plots.display(p)
end


function mesh_test2()
    function p_fn(v)
        return 1.0
    end

    function q_fn(v) 
        k1 = (v.x * 2 ./sqrt(3.) + v.y / sqrt(3)) 
        k2 = v.y 
        return  sqrt(1 + 4* cos(π*sqrt(3) * k1) * cos(π * k2) + 4* cos(π* k2)^2)
    end

    function f_fn(v)
        return q_fn(v)
    end

    t0  = time()
    mesh  = get_square_grid(1., 1., 31)
    surf = Surface2D(mesh, q_fn)
    # println("time: ", time() - t0, " create mesh")
    p = plot(surf)

    cmax = max(maximum(surf._values), 0)
    cmin = min(minimum(surf._values), 0)

    if cmax > 0
        for lnc = -6.:1.: log(cmax)+1.
            mesh  = get_square_grid(1., 1., 57)
            surf = Surface2D(mesh, q_fn)
            level_set_edges = level_set!(surf, exp(lnc), 4, 1e-4)
            plot_level_set!(level_set_edges)
            println("time: ", time() - t0," lnc: ", lnc, " level set: ", length(level_set_edges))
        end
    end


    if cmin < 0 
        for lnc = -6.:1.: log(-cmin)+1.
            mesh  = get_square_grid(1., 1., 57)
            surf = Surface2D(mesh, q_fn)
            level_set_edges = level_set!(surf, -exp(lnc), 4, 1e-4)
            plot_level_set!(level_set_edges)
            println("time: ", time() - t0," lnc: ", lnc, " level set: ", length(level_set_edges))
        end
    end
    Plots.display(p)
end

    

function green_test()

    function ε_k(k)
        return  cos.(2* π *  k.x) + cos.(2 * π * k.y)
    end

    t0  = time()
    Ns = []

    ωs = range(-3., 3.,30)

    for ω in ωs
        BZ = get_square_grid(1., 1., 51)
        dos = calculate_dos(ω, ε_k, BZ)
        push!(Ns, dos)
        println("time: ", time() - t0, " dos: ", dos)
    end

    println("time: ", time() - t0)
    Plots.plot(ωs, Ns, seriestype=:scatter)
    Plots.ylabel!("G(ω)")
end


function green_test_eta()

    function ε_k(k)
        return  cos.(2* π *  k.x) + cos.(2 * π * k.y)
    end

    t0  = time()
    Ns = []

    ωs = range(0., 3.,30)
    η = 1e-5

    for ω in ωs
        BZ = get_square_grid(1., 1., 51)
        dos = η * calculate_dos(ω, k -> ε_k(k)^2 + η^2, BZ)
        push!(Ns, dos)
        println("time: ", time() - t0, " dos: ", dos)
    end

    println("time: ", time() - t0)
    Plots.plot(ωs, Ns, seriestype=:scatter)
    Plots.ylabel!("G(ω)")
end


function dos_graphene_exact(E)
    F = 1. /16. * (abs(E) + 1)^3 * (3 - abs(E))
    K(m) = ellipk(m)

    if abs(E) <= 1
        return abs(E)/π^2 * 1/sqrt(F) * K(abs(E)/F)
    elseif abs(E) <= 3
        return abs(E)/π^2 * 1/sqrt(abs(E)) * K(F/abs(E))
    else 
        return 0.0
    end
end



function green_test2()

    function ε_k(v)
        k1 = (v.x * 2 ./sqrt(3.) + v.y / sqrt(3)) 
        k2 = v.y 
        return  sqrt(1 + 4* cos(π*sqrt(3) * k1) * cos(π * k2) + 4* cos(π* k2)^2)
    end

    t0  = time()
    Ns = []
    Ns_exact = []
    Ns_quad =[]

    ωs = range(-5., 5.,100)

    for ω in ωs
        BZ = get_square_grid(1., 1., 51)
        dos = calculate_dos(ω, ε_k, BZ)
        dos += calculate_dos(ω, k -> -ε_k(k), BZ)
        dos_exact = dos_graphene_exact(ω)
        dos_quad = calculate_dos_quad(ω, ε_k; 
                                        η = 0.0001, 
                                        atol = 1e-10, 
                                        maxevals = 10000000)
        dos_quad += calculate_dos_quad(ω, k -> -ε_k(k); 
                                        η = 0.0001, 
                                        atol = 1e-10, 
                                        maxevals = 10000000)
        println(dos_exact)
        push!(Ns, dos)
        push!(Ns_exact, dos_exact)
        push!(Ns_quad, dos_quad)
        println("time: ", time() - t0, " dos: ", dos)
    end

    println("time: ", time() - t0)
    Plots.plot(ωs, [Ns_exact], label = "exact", color = :black)
    Plots.plot!(ωs, [Ns_quad], label = "quad", seriestype = :scatter, markersize = 5)
    Plots.plot!(ωs, [Ns], label = "level set", seriestype = :scatter, color = :blue4, markersize = 2)
    # Plots.plot(ωs, [log.(abs.(Ns_exact .- Ns)), log.(abs.(Ns_exact .- Ns_quad))], label = ["Level set" "Quad"])
    Plots.xlabel!("ω")
    Plots.ylabel!("Im G(ω)")
end

"""
Baseline
"""

function calculate_dos_quad(ω, ε_k; η = 1e-5, atol = 1e-6, maxevals = 10000)
    integral, error = dblquad((kx, ky) -> η/(η^2 + (ω - ε_k(Vertex2D(kx, ky)))^2.), 0., 1., 0., 1., atol=atol, maxevals = maxevals)
    return integral * 1/π
end


function compare_green()
    function ε_k(v)
        k1 = (v.x * 2 ./sqrt(3.) + v.y / sqrt(3)) 
        k2 = v.y 
        return  sqrt(1 + 4* cos(π*sqrt(3) * k1) * cos(π * k2) + 4* cos(π* k2)^2)
    end

    ω = 2.

    BZ = get_square_grid(1., 1., 51)
    t0 = time()
    dos_LS = calculate_dos(ω, ε_k, BZ; depth_LS  = 10,
                                       depth_int = 10,
                                       tol_LS = 1e-10, 
                                       tol_int = 1e-10)
    println("Level Set time: ", time() - t0, " dos: ", dos_LS)
    t0 = time()
    dos_quad = calculate_dos_quad(ω, ε_k; 
                                    η = 0.00001, 
                                    atol = 1e-10, 
                                    maxevals = 10000000)

    println("quad time: ", time() - t0, " dos: ", dos_quad)
    println("differences: ", abs(dos_LS - dos_quad))
end

function test_runtime_quad(η)
    function ε_k(v)
        k1 = (v.x * 2 ./sqrt(3.) + v.y / sqrt(3)) 
        k2 = v.y 
        return  sqrt(1 + 4* cos(π*sqrt(3) * k1) * cos(π * k2) + 4* cos(π* k2)^2)
    end

    ω = 2.

    BZ = get_square_grid(1., 1., 51)
    dos_exact = dos_graphene_exact(ω)
    errs = []
    times = []

    for maxevals in [10000, 100000, 1000000, 10000000, 100000000]
        println("maxevals", maxevals)
        dos_LS = calculate_dos_quad(ω, ε_k; 
                                    η = η, 
                                    atol = 1e-10, 
                                    maxevals = maxevals)
 
                          
        runtime =  (@belapsed calculate_dos_quad($ω, $ε_k; 
                                                    η = $η, 
                                                    atol = 1e-10, 
                                                    maxevals = $maxevals))
        err = abs(dos_exact - dos_LS)
        push!(times, runtime)
        push!(errs, err)
    end
    Plots.plot!(times, errs ,xaxis=:log, yaxis=:log, label =  "quad η = $η")
    Plots.xaxis!("runtime(s)")
    Plots.yaxis!("error")
end


function test_runtime_LS()
    function ε_k(v)
        k1 = (v.x * 2 ./sqrt(3.) + v.y / sqrt(3)) 
        k2 = v.y 
        return  sqrt(1 + 4* cos(π*sqrt(3) * k1) * cos(π * k2) + 4* cos(π* k2)^2)
    end

    ω = 2.

    BZ = get_square_grid(1., 1., 51)
    dos_exact = dos_graphene_exact(ω)
    errs = []
    times = []
    for depth = 0:12
        println("depth", depth)
        dos_LS = calculate_dos(ω, ε_k, BZ; depth_LS  = depth,
        depth_int = depth,
        tol_LS = 1e-10, 
        tol_int = 1e-10)
                          
        runtime =  (@belapsed calculate_dos($ω, $ε_k, $BZ; depth_LS  = $depth,
                                            depth_int = $depth,
                                            tol_LS = 1e-10, 
                                            tol_int = 1e-10))
        err = abs(dos_exact - dos_LS)
        push!(times, runtime)
        push!(errs, err)
    end
    p = Plots.plot(times, errs ,xaxis=:log, yaxis=:log, label = "level set", color = :blue4)
    Plots.xaxis!("runtime(s)")
    Plots.yaxis!("error")
    return p
end


function green_test_complex()
    function ε_k(k)
        return  cos.(2* π * k.x) + cos.(2 * π * k.y) 
    end

    function v_k(k)
        return  Vertex2DComplex(- 2π * sin(2* π *  k.x), - 2π * sin(2* π *  k.y))
    end
    p = Plots.plot()
    ωs =  -3.001:0.3:3.

    for eta in [0., 0.5, 1.0, 2.0]
        res = []
        for ω in ωs
            println("eta: ", eta, " omega: ", ω)
            mesh  = get_square_grid(1., 1., 57)
            surf = Surface2D(mesh, k -> (ε_k(k) - ω))
            level_set_edges = level_set!(surf, 0.0, 10, 1e-4)
            η(k) = eta 
            depth = 3
            tol = 1e-5

            acc = 0
            for e in level_set_edges
                p1, p2 = e.p1, e.p2
                p1_new = find_zeros(p1, k ->  ε_k(k) - ω, v_k, η, 30)
                p2_new = find_zeros(p2, k ->  ε_k(k) - ω, v_k, η, 30)
                edge_new = Edge(p1_new, p2_new)
                acc += integrate1D(edge_new, k-> 1. /sqrt((v_k(k).x^2 + v_k(k).y^2)), depth, tol)
            end

            push!(res, real(acc))
        end
        Plots.plot!(ωs, res, label = string(eta), seriestype = :scatter)
    end

    display(p)

end

#   TODO: edge case where there is no level set in real 
#   plane but there are some in the complex plane
function green_test_complex2()
    function ε_k(v)
        k1 = (v.x * 2 ./sqrt(3.) + v.y / sqrt(3)) 
        k2 = v.y 
        return  sqrt(1 + 4* cos(π*sqrt(3) * k1) * cos(π * k2) + 4* cos(π* k2)^2)
    end

    function v_k(k)
        k1 = (k.x * 2 ./sqrt(3.) + k.y / sqrt(3)) 
        k2 = k.y 
        ek = ε_k(k)
        return  1/(2. * ek) * Vertex2DComplex(- 4π*sqrt(3)* sin(π*sqrt(3) * k1) * cos(π * k2), 
                                              - π * sin(π * k2) *4* cos(π*sqrt(3) * k1) - 8. * π * cos(π* k2)* sin(π * k2))
    end
    p = Plots.plot()
    ωs =  -4.05:0.1:4.

    for eta in [0., 0.1, 0.2]
        res = []
       
        for ω in ωs
            acc = 0
            # println(ω_)
            for s in [-1. + 0.0im, 1. + 0.0im]
                # println(ω)
                
                mesh  = get_square_grid(1., 1., 57)
                surf = Surface2D(mesh, k -> (s * ε_k(k) - ω))
                level_set_edges = level_set!(surf, 0.0, 7, 1e-10)
                η(k) = eta 
                depth = 5
                tol = 1e-6

                for e in level_set_edges
                    p1, p2 = e.p1, e.p2
                    p1_new = find_zeros(p1, k ->  s * ε_k(k) - ω, k -> s * v_k(k), η, 20)
                    p2_new = find_zeros(p2, k ->  s * ε_k(k) - ω, k -> s * v_k(k), η, 20)
                    edge_new = Edge(p1_new, p2_new)
                    acc += integrate1D(edge_new, k-> 1. /sqrt((v_k(k).x^2 + v_k(k).y^2)), depth, tol)
                end
            end
            println("eta: ", eta, " omega: ", ω, " res: ", acc)
            push!(res, real(acc))
        end
        Plots.plot!(ωs, res, label = string(eta), marker=(:circle,3))
    end

    display(p)

end

function test_runtime()
    p = test_runtime_LS()
    for η in [1e-5, 1e-4, 1e-3, 1e-2]
        println(η)
        test_runtime_quad(η)
    end

    display(p)
end

# default(palette = palette(:PiYG_10))
default(palette = palette(:tab10))
# mesh_test1()
# green_test()
# green_test_eta()
green_test_complex2()



