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

# ╔═╡ 964756ea-6b75-11ef-3685-3fc0ba257201
begin
	import Pkg
	Pkg.activate(temp=true)
	using FEM, Plots, PlutoUI, Roots
end

# ╔═╡ 240fc20e-c154-4188-ac49-e12588737272
using SymPy

# ╔═╡ 181f96f1-e89a-445c-8417-b4ea2b24113f
struct TEDMT{T<:Dim3} <: AbstractConstitutiveLaw{T}
    T0::Float64
    Chom::Function
    κhom::Function
    Fd::Function
    Rd::AbstractRd
    tol::Float64
    function TEDMT{T}(E::Float64, ν::Float64, α::Float64, T0::Real, Rd::AbstractRd; tol::Float64 = 1e-6) where T
        k3 = E / (1 - 2ν)
        μ2 = E / (1 + ν)
        κ = k3 * α
        η1 = 16(1-ν^2) / 9(1-2ν)
        η2 = 32(1-ν) * (5-ν) / 45(2-ν)
        Chom(d) = k3 / (1 + η1*d) * Mandel.J + μ2  / (1 + η2 * d) * Mandel.K
        κhom(d) = κ / (1 + η1*d) * Mandel.δ
        Ta(t) = (t - T0) / log(t / T0)
        Sd = η1 / k3 * Mandel.J + η2 / μ2 * Mandel.K
        Fd(σ, d, t) = begin
            ΔT = t - T0
            σ' * Sd * σ / 2 - 3η1 * α * κ * ΔT / (1 + η1 * d)^2 * (ΔT/2 + Ta(t) - t)
        end
        return new{T}(
            T0, Chom, κhom, Fd, Rd, tol
        )
    end
end

# ╔═╡ 69203f4e-a320-4baa-aefb-fe4caae61e5c
function FEM.constitutive_law_apply!(F::TEDMT{<:Dim3}, p::AbstractPoint{<:AbstractCellType, <:Dim3})
    d, T = p.statev[1:2]
    ΔT = T - F.T0
    ε = p.ε + p.dε
    σ(d) = F.Chom(d) * ε - F.κhom(d) * ΔT
    g(x) = F.Fd(σ(d+x), d+x, T) - Rd(F.Rd, d+x)
    if g(0.) > F.tol
        d += Roots.find_zero(g, 0.0)
    end
    p.statev[1] = d
    p.D .= F.Chom(d)
    p.σ .= σ(d)
end

# ╔═╡ 0f55916b-20fa-4b93-8f1a-fc2fe660428a
function build_con(material, nincr)
	E, ν, α, rc, dc, T0, T1 = material
    ΔT = (T1 - T0) / nincr
	p1 = simple_point_gen(nstatev = 2)
	p1.statev[2] = T0
	constitutive_law_apply!(constitutive_linear_elastic{Dim3}(E, ν), p1)
	cl = TEDMT{Dim3}(E, ν, α, T0, Rd_gen(1, rc, dc))
	rcd = testcon()
	tc = simple_analysis([], [], [], 6, tolerance=1e-5)
	return rcd, tc, p1, cl, ΔT
end

# ╔═╡ d41f5a9a-e5f8-4726-8229-3300264a3555
function thermal_test!(rc::AbstractRecorder, tc::simple_analysis, p::AbstractPoint, cl::TEDMT, nincr, ΔT)
    for incr in 1:nincr
        incr_time = @elapsed begin
            p.statev[2] += ΔT # temperature load
            dε = displacement_load(tc, K = p.D)
            de = dε
            push!(rc.iteration, 0)
            for iter in 1:30
                p.dε .= dε
                constitutive_law_apply!(cl, p)
                tc.f_int .= p.σ
                tc.R .= tc.f_ext - tc.f_int[tc.free_dofs]
                if sqrt( (de'*de) / (dε'*dε) ) <= tc.tolerance
                    rc.iteration[incr] = iter
                    break
                end
                de = displacement_converge(tc, K = p.D)
                dε += de
                rc.iteration[incr] = iter
            end
            p.ε .+= dε
        end
        record!(rc, p)
        push!(rc.time, incr_time)
    end
end

# ╔═╡ 20b9ca20-08aa-416c-962e-60cb72807639
@bind material PlutoUI.combine() do Child
	md"""
	# 材料常数
	E = $(
		Child(NumberField(0:1e6, default=7e4))
	) Mpa, 
	ν = $(
		Child(NumberField(0:0.01:0.5, default=0.2))
	), 
	α = $(
		Child(NumberField(0:1e-5:100, default=2e-5))
	)
	# 抗力函数
	rc = $(
		Child(NumberField(0:1e-5:100, default=9.26e-3))
	), 
	dc = $(
		Child(NumberField(0:1e-5:100, default=7.))
	)
	# 温度加载
	initial temperature = $(
		Child(NumberField(0:1000, default=20))
	)， 
	final temperature = $(
		Child(NumberField(-1000:1000, default=200))
	)
	"""
end

# ╔═╡ 9febd00c-55f4-49da-9025-1a7d67fe2d97
md"""
加载步 = $(@bind nincr Slider(1:1000, show_value=true, default=100))
"""

# ╔═╡ b83efec3-0564-4a7a-b3e5-6e1403d55653
rcd, tc, p1, cl, ΔT = build_con(material, nincr)

# ╔═╡ 1d0ca23e-a080-47f4-a6f2-046ca5df3961
let
	try thermal_test!(rcd, tc, p1, cl, nincr, ΔT) 
	catch
		@info "not converged"
		d = [v[1] for v in rcd.statev]
		T = [v[2] for v in rcd.statev]
		plot(T, d)
	end
	d = [v[1] for v in rcd.statev]
	T = [v[2] for v in rcd.statev]
	plot(T, d, xlabel="T", ylabel="d", title="T-d curve")
end

# ╔═╡ 05df9daf-0762-4e22-9179-39567e7ae1de
rcd

# ╔═╡ 7a481c86-6dfc-46fc-a40c-aab6f08115de
@syms t, t0

# ╔═╡ fb4500b1-6b56-41d1-afe7-09fb3910338a
δt = t - t0

# ╔═╡ a5627de6-f6f1-4e1b-b0b3-551dfa1b6e99
ta = δt / log(t/t0)

# ╔═╡ bd31105e-b175-44a6-8f36-454836b35a4d
y = δt(δt / 2 + ta - t)

# ╔═╡ 6a063df5-4200-495d-9cf4-0c6d04fe9639
simplify(y)

# ╔═╡ 94af2024-eccf-4056-9750-f521abc378ac
z = y.subs(t0, 200)

# ╔═╡ 08abb2c1-f0db-4f67-a4c3-e44cdd73ef10
t1 = collect(0:1000)

# ╔═╡ 34c048e3-1208-435e-8e69-85a06033a1b1
z1 = z.subs.(t, t1)

# ╔═╡ aa2eb9ae-c6db-4d0f-a97e-2585b978743c
plot(t1, z1)

# ╔═╡ 2eb8bd26-a842-4395-820f-19514ef8be79


# ╔═╡ Cell order:
# ╠═964756ea-6b75-11ef-3685-3fc0ba257201
# ╠═181f96f1-e89a-445c-8417-b4ea2b24113f
# ╠═69203f4e-a320-4baa-aefb-fe4caae61e5c
# ╠═0f55916b-20fa-4b93-8f1a-fc2fe660428a
# ╠═b83efec3-0564-4a7a-b3e5-6e1403d55653
# ╠═d41f5a9a-e5f8-4726-8229-3300264a3555
# ╟─20b9ca20-08aa-416c-962e-60cb72807639
# ╟─9febd00c-55f4-49da-9025-1a7d67fe2d97
# ╠═1d0ca23e-a080-47f4-a6f2-046ca5df3961
# ╠═05df9daf-0762-4e22-9179-39567e7ae1de
# ╠═240fc20e-c154-4188-ac49-e12588737272
# ╠═7a481c86-6dfc-46fc-a40c-aab6f08115de
# ╠═fb4500b1-6b56-41d1-afe7-09fb3910338a
# ╠═a5627de6-f6f1-4e1b-b0b3-551dfa1b6e99
# ╠═bd31105e-b175-44a6-8f36-454836b35a4d
# ╠═6a063df5-4200-495d-9cf4-0c6d04fe9639
# ╠═94af2024-eccf-4056-9750-f521abc378ac
# ╠═08abb2c1-f0db-4f67-a4c3-e44cdd73ef10
# ╠═34c048e3-1208-435e-8e69-85a06033a1b1
# ╠═aa2eb9ae-c6db-4d0f-a97e-2585b978743c
# ╠═2eb8bd26-a842-4395-820f-19514ef8be79
