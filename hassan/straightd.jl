using FEM, DelimitedFiles
using Debugger

struct straightd{T<:Dim3} <: AbstractConstitutiveLaw{T}
    a::Float64
    x::Int
    y::Int
    k3::Float64
    μ2::Float64
    weaken_ratio::Dict{Bool, NTuple{4, Function}}
    isopen::Function
    Rd::AbstractRd
    tol::Float64

    function straightd{T}(S::Symbol, E::Float64, ν::Float64, Rd::AbstractRd; tol::Float64 = 1e-6, l::Float64) where T
        @assert S == :MT "only MT"
        m = elastic_damage{T}(S, E, ν, Rd, tol= tol, l = l)
        return new{T}(m.a, 1, 1, m.k3, m.μ2, m.weaken_ratio, m.isopen, m.Rd, m.tol)
    end
end
function FEM.constitutive_law_apply!(F::straightd{T}, p::AbstractPoint{<:AbstractCellType, T}) where T<:Dim3
    d = p.statev[1]
    ε = p.ε + p.dε
    kd, dkd, μd, dμd = F.weaken_ratio[F.isopen(p.D * ε)]
    k3 = F.k3
    μ2 = F.μ2
    g(x) = -ε' * (dkd(d+x) * k3 * Mandel.J + dμd(d+x) * μ2 * Mandel.K) * ε / 2 - Rd(F.Rd, d+x)
    if g(0.) > F.tol
        d += Roots.find_zero(g, 0.0)
    end
    Chom = kd(d) * k3 * Mandel.J + μd(d) * μ2 * Mandel.K
    p.D .= Chom
    p.σ .= Chom * ε
    p.statev[1] = d
    return nothing
end

function model_gen(nincr)
    e1 = readdlm("elements1", Int)
    e2 = readdlm("elements2", Int)
    e3 = readdlm("elements3", Int)
    e = hcat(e1', e2', e3')
    nodes = readdlm("nodes")' # 单位mm

    E = 7e4
    ν = 0.15
    rc = 9.36e-3
    dc = 7.0
    l = 3.0 # length scale 单位mm
    cl = straightd{Dim3}(:MT, E, ν, Rd_gen(1, rc, dc), l = 3.)

    mesh1 = fem_mesh_gen(nodes)
    fem1 = elastoplastic(H1, e, nodes, 3, 3, 3, E = E, ν = ν, nstatev = 3)
    rcd = fem_recorder("st1", "output_g", fem1.quadrature_points)

    top = findall(nodes[2, :].==maximum(nodes[2, :]))
    bottom = findall(nodes[2, :].==minimum(nodes[2, :]))
    fixed_dofs = vcat(3bottom.-2, 3bottom.-1, 3bottom) # 底部固定
    load_dofs = 3top.-1 # 顶部y方向加载
    max_load = 4e-1
    load_per_step = max_load / nincr
    load_value = fill(load_per_step, length(load_dofs))

    solve1 = simple_analysis(fixed_dofs, load_dofs, load_value, mesh1.dofs, tolerance = 1e-6)
    return rcd, cl, mesh1, fem1, solve1
end

function analysis!(rc::AbstractRecorder,
                   a::simple_analysis,
                   fe::AbstractFem,
                   mesh::AbstractMesh,
                   g::AbstractConstitutiveLaw;
                   nincr::Int,
                   miter::Int = 30,
                   u::Vector{Float64} = zeros(mesh.dofs),
                   K_u = elastic_initialization!(fe, mesh),
                   K_g = gradient_stiffness_assemble(g, fe, mesh)
                   )

    for incr in 1:nincr
        incr_time = @elapsed begin
            du = displacement_load(a, K = K_u)
            converged = false
            for iter in 1:miter
                displacement_apply!(fe, du = du)
                constitutive_law_apply!(g, fe)
                rhs = K_g \ gradient_rhs_gen(g, fe, mesh)
                gradient_apply!(g, fe, rhs)
                K_u, f_int = elastic_system_assemble!(fe, mesh)
                u += du
                displacement_apply!(fe, u = u)
                a.f_int .= f_int
                a.R .= a.f_ext - a.f_int[a.free_dofs]
                du = displacement_converge(a, K = K_u)
                if norm(du) / norm(u) <= a.tolerance && norm(a.R) / norm(a.f_int) <= a.tolerance
                    converged = true
                    push!(rc.iteration, iter)
                    print("Incriment $incr converged at $iter", " th iteration. ")
                    break
                end
            end
            if !converged
                push!(rc.iteration, miter)
                @warn "The solution did not converge within $miter iterations."
            end
            push!(rc.convergence, converged)
            record!(rc, a, fe, u, i = incr)
        end
        push!(rc.time, incr_time)
        println("Time comsume: $incr_time", "s.")
    end
end



function analysis(nincr)
    rcd, cl, mesh1, fem1, solve1 = model_gen(nincr)
    analysis!(rcd, solve1, fem1, mesh1, cl, nincr = nincr, miter = 100)
    record!(rcd)
    return rcd, mesh1, fem1
end

function output(nincr, rcd, mesh, fe)
    for i in 1:nincr
        plot2vtu(rcd, mesh, i, fe.element_table)
        println(i, " done.")
    end
end
