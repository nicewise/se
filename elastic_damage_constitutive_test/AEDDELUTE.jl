### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
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
using StaticArrays

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
struct AEDDELUTE{T<:Dim3} <: AbstractConstitutiveLaw{T}
    C0::SMatrix{6, 6, Float64, 36}
    cn::Float64
    ct::Float64
    n::Int # 代表性裂隙族数量，ie. 球面积分高斯点数目
    cracks::Vector{MicroCrack}
    Chom::Function
    isopen::Function
    Rd::AbstractRd
    tol::Float64
    function AEDDELUTE{T}(E::Float64, ν::Float64, Rd::AbstractRd; tol::Float64 = 1e-6, n::Int) where T
        known_quadrature = [21, 33, 48, 96, 99, 144, 198]
        @assert n ∈ known_quadrature "只知道($i for i in known_quadrature), 不知道$n 积分点怎么办."
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
        κe = Dict(
            true => d -> d < 3(1-2ν)/16(1-ν)^2 ?
                (-16*d*(ν - 1)^2 - 6*ν + 3)/(32*d*ν^2*(ν - 1) - 6*ν + 3) : 0,
            false => d -> 1
        )
        κμ(d) = d < 3(2-ν)/16(1-ν) ? 
            1 - 16(1-ν) / 3(2-ν) * d : 0
        Chom(d::Real, key::Bool) = begin
            κ = κe[key]
            r = 1 - ν - 2κ(d) * ν^2
            return [
                E / r
                E * κ(d) * (1-ν) / r
                E / (1 + ν)
                2κμ(d)
                E * κ(d) * ν / r
                E * κ(d) * ν / r
            ]
        end
        isopen = σ -> [σ' * crack.N for crack in cracks] .>= 0 # σ -> BitVector
        return new{T}(C0, cn, ct, n, cracks, Chom, isopen, Rd, tol)
    end
end

# ╔═╡ 46616fdb-de38-4116-9a28-cfc0fb003f65
Fd(F::AEDDELUTE{Dim3}, ε::Vector{Float64}, ckeys::BitVector) =
    # keys to determine the specific micro crack is open or not. ckeys = crack_keys
    [- ε' * F.C0 * Sn(F.cracks[i], key, F.cn, F.ct) * F.C0 * ε / 2 for (i, key) in enumerate(ckeys)]

# ╔═╡ b2fc2e95-9b0e-4c9a-9241-5bcf83bf69f4
function FEM.constitutive_law_apply!(F::AEDDELUTE{T}, p::AbstractPoint{<:AbstractCellType, T}) where T<:Dim3
    ε = p.ε + p.dε
    d = p.statev[1:F.n]
    stats = F.isopen(p.D * ε) # Vector{Bool} indicate is microcrack open or not
    Fd_const = Fd(F, ε, stats)
    g = [x->Fd_const[i] - Rd(F.Rd, d[i]+x) for i in eachindex(d)]
    #g = [x->Fd_const[i] - F.cracks[i].ω * Rd(F.Rd, d[i]+x) for i in eachindex(d)]
    for (i, f) in enumerate(g)
        if f(0.) > F.tol
            d[i] += Roots.find_zero(f, 0.0)
        end
    end
    p.statev[1:F.n] .= d
    Chom = sum(compute_walpole(F.cracks[i], F.Chom(d[i], stats[i])) * F.cracks[i].ω for i in 1:F.n)
    p.D .= Chom
    p.σ .= Chom * ε
    return nothing
end

# ╔═╡ b9e69b5f-ae1e-42ee-9402-3c8e15dd5316
E, ν = 7e4, 0.2

# ╔═╡ 9bc569db-9b73-4d55-aedf-6c04d6f5c4b7
function build_con(material, max_load, nincr)
	E, ν, rc, dc = material
	ng = 33
	load_per_incr = max_load / nincr
	p = simple_point_gen(nstatev = ng)
	constitutive_law_apply!(constitutive_linear_elastic{Dim3}(E, ν), p)
	cl = AEDDELUTE{Dim3}(E, ν, Rd_gen(1, rc, dc), n = ng)
	rcd = testcon()
	tc = simple_analysis([], [1], [load_per_incr], 6, tolerance=1e-5)
	return rcd, tc, p, cl
end

# ╔═╡ 92f9d2a9-465c-49b4-915e-65ed8f177086
@bind material PlutoUI.combine() do Child
	md"""
	# 材料常数
	E = $(
		Child(NumberField(0:1e6, default=7e4))
	) Mpa, 
	ν = $(
		Child(NumberField(0:0.01:0.5, default=0.2))
	), 

	# 抗力函数
	rc = $(
		Child(NumberField(0:1e-5:100, default=9.26e-3))
	), 
	dc = $(
		Child(NumberField(0:1e-5:100, default=7.))
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
	# plot(ε1, d)
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
# ╠═46616fdb-de38-4116-9a28-cfc0fb003f65
# ╠═b2fc2e95-9b0e-4c9a-9241-5bcf83bf69f4
# ╠═b9e69b5f-ae1e-42ee-9402-3c8e15dd5316
# ╠═9bc569db-9b73-4d55-aedf-6c04d6f5c4b7
# ╠═fb412e76-34a5-4388-821a-80af1f89c2e5
# ╟─92f9d2a9-465c-49b4-915e-65ed8f177086
# ╟─08997cb0-0ed4-4f84-89f3-72bc396ee210
# ╟─e38f14ed-153a-4543-839b-8ae84c9f8739
# ╠═4a32d743-c37a-407b-8b11-23fd5de71e77
# ╠═3f5a394c-4a08-47d1-9d2e-8a9ba63263a5
