### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ c400698c-6a9c-11ef-22c3-21215802f184
begin
	import Pkg
	Pkg.activate(temp=true)
	using FEM, Plots, PlutoUI, Roots
end

# ╔═╡ b4da2b00-1a3a-48f4-b419-adc0cc2efa6c
using StaticArrays, LinearAlgebra

# ╔═╡ 42f2bdfc-f3e9-4be6-9dbd-bf06c0919e0d
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

# ╔═╡ 8db8a60d-a30d-427e-a868-6b6947cb0d82
compute_walpole(c::MicroCrack, x) =
    x[1] * c.E1 + x[2] * c.E2 + x[3] * c.E3 + x[4] * c.E4 + x[5] * c.E5 + x[6] * c.E6

# ╔═╡ db49867e-cd0e-4606-9c71-f9a75f0a163f
@inline function Sn(crack::MicroCrack, key::Bool, cn, ct) # cn和ct是书P60页所述cn、ct之倒数，为了少做除法
    if key
        return cn * crack.E2 + ct * crack.E4
    else
        return ct * crack.E4
    end
end

# ╔═╡ 3e9a5fe8-ec9e-4695-a51d-decbe1d13ed1
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
        return new{T}(l, x, y, n, dnw[:, 4], cracks, C0, Tr_gen, Pd, isopen, Rd, tol)
    end
end

# ╔═╡ 1f80e444-54b5-4904-aba0-2bf7117f7193
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

# ╔═╡ b2fc2e95-9b0e-4c9a-9241-5bcf83bf69f4
function FEM.constitutive_law_apply!(F::AEDPCW{T}, p::AbstractPoint{<:AbstractCellType, T}; before_gradient::Union{Nothing, Bool}=nothing) where T<:Dim3
    if isnothing(before_gradient) || before_gradient
        d_old = p.statev[F.x]
        ε = p.ε + p.dε
        d(x) = d_old + x # vector -> vector
        stats = F.isopen(p.D * ε) # Vector{Bool} indicate is microcrack open or not
        Tr = [F.Tr_gen(crack, stats[i]) for (i, crack) in enumerate(F.cracks)]
        Cd(x) = sum(d(x) .* Tr) # vector -> matrix
        B(x) = (Mandel.I + F.Pd * Cd(x)) \ F.Pd # vector -> matrix
        CdB(x) = 2 * B(x) * Cd(x) - Mandel.I # vector -> matrix
        Fd(x) = [-ε'*Tr[i]*CdB(x)*ε/2 for i in 1:F.n]
        #g(x) = Fd(x) - Rd(F.Rd, d(x)) # vector -> vector
        g(x) = Rd(F.Rd, d(x), F.ω) - Fd(x)# vector -> vector
        if sum( g(zeros(F.n)) .< -F.tol ) >= 1
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
        Tr = [F.Tr_gen(crack, stats[i]) for (i, crack) in enumerate(F.cracks)]
        cd = sum(d .* Tr)
        b = (Mandel.I + F.Pd * cd) \ F.Pd
        Chom = F.C0 - cd + cd * b * cd
        p.D .= Chom
        p.σ .= Chom * ε
    end
    return nothing
end

# ╔═╡ 9bc569db-9b73-4d55-aedf-6c04d6f5c4b7
function build_con(material, max_load, nincr)
	E, ν, c0, c1 = material
	ng = 33
	load_per_incr = max_load / nincr
	p = simple_point_gen(nstatev = ng)
	constitutive_law_apply!(constitutive_linear_elastic{Dim3}(E, ν), p)
	cl = AEDPCW{Dim3}(E, ν, Rd_gen(4, c0, c1), n = ng)
	rcd = testcon()
	tc = simple_analysis([], [1], [load_per_incr], 6, tolerance=1e-5)
	return rcd, tc, p, cl
end

# ╔═╡ 0efa2aee-27c6-4fb9-bf0c-f808ab3087e0
Rd_gen(4, 7.4e-4, 4e-3)

# ╔═╡ 92f9d2a9-465c-49b4-915e-65ed8f177086
@bind material PlutoUI.combine() do Child
	md"""
	# 材料常数
	E = $(
		Child(NumberField(0:1e6, default=3.8e4))
	) Mpa, 
	ν = $(
		Child(NumberField(0:0.01:0.5, default=0.19))
	), 

	# 抗力函数
	c0 = $(
		Child(NumberField(0:1e-5:100, default=7.4e-4))
	), 
	c1 = $(
		Child(NumberField(0:1e-5:100, default=4e-3))
	)
	"""
end

# ╔═╡ 08997cb0-0ed4-4f84-89f3-72bc396ee210
md"""
应变 = $(@bind max_load Slider(-2e-2:1e-3:2e-1, show_value=true, default=2e-2))
"""

# ╔═╡ e38f14ed-153a-4543-839b-8ae84c9f8739
md"""
加载步 = $(@bind nincr Slider(1:1000, show_value=true, default=100))
"""

# ╔═╡ fb412e76-34a5-4388-821a-80af1f89c2e5
rcd, tc, p, cl = build_con(material, max_load, nincr)

# ╔═╡ 4a32d743-c37a-407b-8b11-23fd5de71e77
let
	testcon!(rcd, tc, p, cl, nincr)

	σ1 = vcat(0, [σ[1] for σ in rcd.σ])
	ε1 = vcat(0, [ε[1] for ε in rcd.ε])
	ε2 = vcat(0, [ε[2] for ε in rcd.ε])
	d = vcat(0, [v[1] for v in rcd.statev])
	plot(ε1, σ1)
	plot!(ε2, σ1)
	#plot(ε1, d)
end

# ╔═╡ 3f5a394c-4a08-47d1-9d2e-8a9ba63263a5
rcd

# ╔═╡ Cell order:
# ╟─c400698c-6a9c-11ef-22c3-21215802f184
# ╠═b4da2b00-1a3a-48f4-b419-adc0cc2efa6c
# ╠═42f2bdfc-f3e9-4be6-9dbd-bf06c0919e0d
# ╠═8db8a60d-a30d-427e-a868-6b6947cb0d82
# ╠═db49867e-cd0e-4606-9c71-f9a75f0a163f
# ╠═3e9a5fe8-ec9e-4695-a51d-decbe1d13ed1
# ╠═1f80e444-54b5-4904-aba0-2bf7117f7193
# ╠═b2fc2e95-9b0e-4c9a-9241-5bcf83bf69f4
# ╠═9bc569db-9b73-4d55-aedf-6c04d6f5c4b7
# ╠═fb412e76-34a5-4388-821a-80af1f89c2e5
# ╠═0efa2aee-27c6-4fb9-bf0c-f808ab3087e0
# ╟─92f9d2a9-465c-49b4-915e-65ed8f177086
# ╟─08997cb0-0ed4-4f84-89f3-72bc396ee210
# ╟─e38f14ed-153a-4543-839b-8ae84c9f8739
# ╠═4a32d743-c37a-407b-8b11-23fd5de71e77
# ╠═3f5a394c-4a08-47d1-9d2e-8a9ba63263a5
