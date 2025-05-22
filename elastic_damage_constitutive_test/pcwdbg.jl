using FEM, DelimitedFiles
using Debugger
using LinearAlgebra
using StaticArrays
using Plots
using JuMP, PATHSolver

struct MicroCrack
    direction::SVector{3}
    ω::Float64 # weight
    N::SVector{6, Float64} # 法向方向张量
    E1::SMatrix{6, 6, Float64} # Walpole 张量基
    E2::SMatrix{6, 6, Float64}
    E3::SMatrix{6, 6, Float64}
    E4::SMatrix{6, 6, Float64}
    E5::SMatrix{6, 6, Float64}
    E6::SMatrix{6, 6, Float64}
    function MicroCrack(n::AbstractVector{<:Real}, w::Float64)
        @assert length(n) == 3 "Input vector 'n' must have length 3, got $(length(n))"
        return MicroCrack(SVector{3}(n), w)
    end
    function MicroCrack(n::SVector{3}, w::Float64)
        N = n ⊗ n
        T = Mandel.δ - N
        return new(
            n,
            w,
            N,
            T ⊗ T / 2,
            N ⊗ N,
            T ⋆ T - T ⊗ T / 2, # ⋆ = ⊗s
            N ⋆ T + T ⋆ N,
            N ⊗ T,
            T ⊗ N
        )
    end
end

compute_walpole(c::MicroCrack, x) =
    x[1] * c.E1 + x[2] * c.E2 + x[3] * c.E3 + x[4] * c.E4 + x[5] * c.E5 + x[6] * c.E6

@inline function Sn(crack::MicroCrack, key::Bool, cn, ct) # cn和ct是书P60页所述cn、ct之倒数，为了少做除法
    if key
        return cn * crack.E2 + ct * crack.E4
    else
        return ct * crack.E4
    end
end

struct LCP
    n::Int
    model
    function LCP(n) # n denotes the number of variables
        model = Model(PATHSolver.Optimizer)
        #set_optimizer_attribute(model, "output", "no")
        set_optimizer_attribute(model, "silent", true)
        # without this option set to "true"
        # PATHSolver will print a lot of information to stdout
        return new(n, model)
    end
end
function LCP_solve(lcp::LCP, M, q)
    x = @variable(lcp.model, [1:lcp.n], lower_bound = 0)
    @constraint(lcp.model, M * x .+ q ⟂ x)
    optimize!(lcp.model)
    return(value.(x)) # obtain solution
end

struct AEDPCW{T<:Dim3} <: AbstractConstitutiveLaw{T}
    l::Float64
    x::Vector{Int}
    y::Vector{Int}
    n::Int # 代表性裂隙族数量，ie. 球面积分高斯点数目
    ω::Vector{Float64}
    cracks::Vector{MicroCrack}
    C0::SMatrix{6, 6, Float64}
    Tr_gen::Function
    Pd::SMatrix{6, 6, Float64}
    isopen::Function
    Rd::AbstractRd
    tol::Float64
    lcp::LCP
    function AEDPCW{T}(E::Float64, ν::Float64, Rd::AbstractRd; tol::Float64 = 1e-6, l::Float64 = 0., n::Int) where T
        x = collect(1:n)
        y = collect(n+1:2n)
        dnw = sphere_quad(n) # direction and weight
        cracks = MicroCrack[]
        for i in axes(dnw, 1)
            push!(cracks, MicroCrack(dnw[i, 1:3], dnw[i, 4]))
        end
        cn = 16(1-ν^2) / 3E
        ct = (2-ν) * cn
        k3 = E / (1 - 2ν)
        μ2 = E / (1 + ν)
        C0 = k3 * Mandel.J + μ2 * Mandel.K
        Tr_gen(crack::MicroCrack, key::Bool) = C0 * Sn(crack, key, cn, ct) * C0 * crack.ω
        αk = 1 / (k3 + 2μ2)
        αμ = (2k3 + 6μ2) / 5(k3 + 2μ2) / μ2
        Pd = αk * Mandel.J + αμ * Mandel.K
        isopen = σ -> [σ' * crack.N for crack in cracks] .>= 0 # σ -> vector{Bool}
        lcp = LCP(n)
        return new{T}(l, x, y, n, dnw[:, 4], cracks, C0, Tr_gen, Pd, isopen, Rd, tol, lcp)
    end
end

function NCP_desent(f, nv::Int)
    α = 0.01
    G = 0.01
    H(x) = max.(x - f(x) / G, 0)
    F(x) = -f(x)' * (H(x)-x) - (H(x)-x)' * (H(x)-x) * G / 2
    z = zeros(nv)
    β = 1
    miter = 100
    tol = 1e-6
    for iter in 1:miter
        Fn = F(z)
        d = H(z) - z
        while F(z+β*d) - F(z) > α * β * norm(z)
            β /= 2
            if β < tol / 2
                #print(β)
                break
            end
        end
        z += β*d
        if norm(F(z) - Fn) / norm(Fn) <= tol
            @info "求解NCP迭代了 $iter 步"
            break
        end
        if iter == miter
            @warn "NCP 没解出来"
        end
    end
    return(z)
end




function FEM.constitutive_law_apply!(F::AEDPCW{T}, p::AbstractPoint{<:AbstractCellType, T}; before_gradient::Union{Nothing, Bool}=nothing) where T<:Dim3
    if isnothing(before_gradient) || before_gradient
        d_old = p.statev[F.x]
        ε = p.ε + p.dε
        d(x) = d_old + x # vector -> vector
        stats = F.isopen(p.D * ε) # Vector{Bool} indicate is microcrack open or not
        Tr = [F.Tr_gen(crack, stats[i]) for (i, crack) in enumerate(F.cracks)]
        Cd(x) = sum(d(x) .* Tr) # vector -> matrix
        B(x) = (Mandel.I + F.Pd * Cd(x)) \ F.Pd # vector -> matrix
        BCd2(x) = 2 * B(x) * Cd(x) - Mandel.I # vector -> matrix
        Fd(x) = [-ε'*Tr[i]*BCd2(x)*ε/2 for i in 1:F.n]
        g(x) = Fd(x) - Rd(F.Rd, d(x), F.ω) # vector -> vector

        CdB(x) = Mandel.I - Cd(x) * B(x)
        BCd(x) = Cd(x) * B(x) - Mandel.I
        function dFd(x)
            o = zeros(F.n, F.n)
            for i in 1:F.n
                o[i, i] = dFd_val(x, Tr[i], Tr[i]) / 2
                for j in i+1:F.n
                    o[i, j] = dFd_val(x, Tr[j], Tr[i])
                end
            end
            return (o + o')
        end
        dFd_val(x, Tj, Ti) = -ε'*CdB(x)*(Tj*B(x)*Ti)*BCd(x)*ε #TODO

        g∂d = dFd(d_old) - dRd(F.Rd, d_old, F.ω)
        g_penetrator = g(zeros(F.n))
        if sum( g_penetrator .> F.tol ) >= 1
            #Δd = NCP_desent(g, F.n)#optim(g) # TODO
            Δd = LCP_solve(F.lcp, g∂d, g_penetrator)
            d_old = d(Δd)
        end
        p.statev[F.x] .= d_old
    end
    if isnothing(before_gradient) || !before_gradient
        d = isnothing(before_gradient) ? p.statev[F.x] : p.statev[F.y]
        d = p.statev[F.x]
        ε = p.ε + p.dε
        stats = F.isopen(p.D * ε) # Vector{Bool} indicate is microcrack open or not
        Tr = [F.Tr_gen(crack, stats[i]) for (i, crack) in enumerate(F.cracks)]
        cd = sum(d .* Tr)
        b = (Mandel.I + F.Pd * cd) \ F.Pd
        Chom = F.C0 - cd + cd * b * cd
        p.D .= Chom
        p.σ .= Chom * ε
    end
    return nothing
end

function build_con(nincr)
    E, ν, c0, c1 = [3.8e4, 0.19, 7.4e-4, 4e-3]
    max_load = 2e-2
    ng = 33
    load_per_incr = max_load / nincr
    p = simple_point_gen(nstatev = ng)
    constitutive_law_apply!(constitutive_linear_elastic{Dim3}(E, ν), p)
    cl = AEDPCW{Dim3}(E, ν, Rd_gen(4, c0, c1), n = ng)
    rcd = testcon()
    tc = simple_analysis([], [1], [load_per_incr], 6, tolerance=1e-5)
    return rcd, tc, p, cl
end

function analysis(nincr)
    rcd, tc, p, cl = build_con(nincr)
    testcon!(rcd, tc, p, cl, nincr)
    σ1 = vcat(0, [σ[1] for σ in rcd.σ])
    ε1 = vcat(0, [ε[1] for ε in rcd.ε])
    ε2 = vcat(0, [ε[2] for ε in rcd.ε])
    d = vcat(0, [v[1] for v in rcd.statev])
    plot(ε1, σ1)
    plot!(ε2, σ1)
    #plot(ε1, d)
end
