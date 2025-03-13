using FEM, DelimitedFiles
using Debugger
using LinearAlgebra
using StaticArrays

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

struct AEDPCW{T<:Dim3} <: AbstractConstitutiveLaw{T}
    l::Float64
    x::Vector{Int}
    y::Vector{Int}
    n::Int # 代表性裂隙族数量，ie. 球面积分高斯点数目
    cracks::Vector{MicroCrack}
    C0::SMatrix{6, 6, Float64}
    Tr::Dict # FIXME
    Pd::SMatrix{6, 6, Float64}
    isopen::Function
    Rd::AbstractRd
    tol::Float64
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
        Tr = Dict(
            true => c -> C0 * (cn * c.E2 + ct * c.E4) * C0 * c.ω,
            false => c -> C0 * ct * c.E4 * C0 * c.ω
        ) # c -- > crack 
        αk = 1 / (k3 + 2μ2)
        αμ = (2k3 + 6μ2) / 5(k3 + 2μ2) / μ2
        Pd = αk * Mandel.J + αμ * Mandel.K
        isopen = σ -> [σ' * crack.N for crack in cracks] .>= 0 # σ -> vector{Bool}
        return new{T}(l, x, y, n, cracks, C0, Tr, Pd, isopen, Rd, tol)
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
        stats = F.isopen(p.D * ε) # Vector{Bool} indicate is microcrack open or not
        Tr = [F.Tr[stats[i]](F.cracks[i]) for i in 1:F.n] # FIXME
        d(x) = d_old + x # vector -> vector
        Cd(x) = sum(d(x) .* Tr) # vector -> matrix
        B(x) = (Mandel.I + F.Pd * Cd(x)) \ F.Pd # vector -> matrix
        #∂Chom(x) = begin # vector -> vecctor{matrix}
        #    cd = Cd(x)
        #    b = B(x)
        #    return 2Tr .* b * cd - cd * b .* Tr .* b * cd
        #end
        #Fd(x) = - ε' .* ∂Chom(x) .* ε / 2 # vector -> vector

        Fd(x) = begin
            cd = Cd(x)
            b = B(x)
            CdB(x) = 2 * B(x) * Cd(x) - Mandel.I
            o = zeros(F.n)
            Threads.@threads for i in eachindex(Tr)
                o[i] = -ε'*Tr[i]*CdB(x)*ε/2
            end
            o
        end

        g(x) = Fd(x) - Rd(F.Rd, d(x)) # vector -> vector
        if sum( g(zeros(F.n)) .> F.tol ) >= 1
            Δd = NCP_desent(g, F.n)#optim(g) # TODO
            d_old = d(Δd)
        end
        p.statev[F.x] .= d_old
    end
    if isnothing(before_gradient) || !before_gradient
        d = isnothing(before_gradient) ? p.statev[F.x] : p.statev[F.y]
        d = p.statev[F.x]
        ε = p.ε + p.dε
        stats = F.isopen(p.D * ε) # Vector{Bool} indicate is microcrack open or not
        Tr = [F.Tr[stats[i]](F.cracks[i]) for i in 1:F.n] # FIXMTE
        cd = sum(d .* Tr)
        b = (Mandel.I + F.Pd * cd) \ F.Pd
        Chom = F.C0 - cd + cd * b * cd
        p.D .= Chom
        p.σ .= Chom * ε
    end
    return nothing
end


function build_con(nincr)
    E, ν, rc, dc = [7e4, 0.2, 9.26e-3, 7.0]
    max_load = 2e-2
    ng = 33
    load_per_incr = max_load / nincr
    p = simple_point_gen(nstatev = ng)
    constitutive_law_apply!(constitutive_linear_elastic{Dim3}(E, ν), p)
    cl = AEDPCW{Dim3}(E, ν, Rd_gen(1, rc, dc), n = ng)
    rcd = testcon()
    tc = simple_analysis([], [1], [load_per_incr], 6, tolerance=1e-5)
    return rcd, tc, p, cl
end

function analysis(nincr)
    rcd, tc, p, cl = build_con(nincr)
    testcon!(rcd, tc, p, cl, nincr)
    rcd
end
