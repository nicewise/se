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

# ╔═╡ 027dad22-6096-11ef-1d9f-af6c9561225a
begin
	import Pkg
	Pkg.activate(temp=true)
	using FEM, Plots, PlutoUI
end

# ╔═╡ bf1fce26-8530-48b5-9cc4-34fae8adffb4
begin
struct R_energy <: AbstractRd
	V::Float64
	γ::Float64
	N::Float64
	R_energy(V, γ, d0, a0) = new(V, γ, d0/a0^3)
end
FEM.Rd(R::R_energy, d) = 4pi/3*R.γ/(d/R.N)^(1/3)*R.V
x = collect(0:0.1:1000)
R1 = R_energy(1., 2.7e-8, 1e-3, 0.01)
plot(x, [Rd(R1, i) for i in x])
end

# ╔═╡ e9963c5c-3598-43b2-a108-219867bd08ec
function build_con(S, max_load, nincr)
	E = 7e4
	ν = 0.15
	rc = 1.
    dc = 0.03
	load_per_incr = max_load / nincr
	p1 = simple_point_gen(nstatev = 3)
	constitutive_law_apply!(constitutive_linear_elastic{Dim3}(E, ν), p1)
	cl = elastic_damage{Dim3}(S, E, ν, Rd_gen(4, rc, dc))
	rcd = testcon()
	tc = simple_analysis([], [1], [load_per_incr], 6, tolerance=1e-5)
	return rcd, tc, p1, cl
end

# ╔═╡ 7a79e662-63ba-4702-b6a0-e487ea90cd87
md"""
S = $(@bind S Slider([:DELUTE, :MT, :PCW], show_value=true, default=:MT))
"""

# ╔═╡ b36ee995-d5cf-4169-80ce-c9b32447ac55
md"""
应变 = $(@bind max_load Slider(-2e-2:1e-3:2e-2, show_value=true, default=2e-2))
"""

# ╔═╡ f350639b-4314-493f-b315-477a4722d9a7
md"""
加载步 = $(@bind nincr Slider(1:1000, show_value=true, default=100))
"""

# ╔═╡ 1811384b-92cf-4d20-83da-f70f24108124
rcd, tc, p1, cl = build_con(S, max_load, nincr)

# ╔═╡ 95bd0f1b-889a-457e-a709-9fe753fda30c
let
	testcon!(rcd, tc, p1, cl, nincr)
	
	σ1 = vcat(0, [σ[1] for σ in rcd.σ])
	ε1 = vcat(0, [ε[1] for ε in rcd.ε])
	ε2 = vcat(0, [ε[2] for ε in rcd.ε])
	d = [v[1] for v in rcd.statev]
	plot(ε1, σ1)
	plot!(ε2, σ1)
	# plot!(d)
end

# ╔═╡ f15ebcd9-0755-4c9f-a417-21eb20e707c8
testcon!(rcd, tc, p1, cl, 100)

# ╔═╡ cc5e59c1-f39f-40b6-b9b8-58d3adb8a54e
rcd

# ╔═╡ Cell order:
# ╟─027dad22-6096-11ef-1d9f-af6c9561225a
# ╠═bf1fce26-8530-48b5-9cc4-34fae8adffb4
# ╠═e9963c5c-3598-43b2-a108-219867bd08ec
# ╠═1811384b-92cf-4d20-83da-f70f24108124
# ╟─7a79e662-63ba-4702-b6a0-e487ea90cd87
# ╟─b36ee995-d5cf-4169-80ce-c9b32447ac55
# ╟─f350639b-4314-493f-b315-477a4722d9a7
# ╠═95bd0f1b-889a-457e-a709-9fe753fda30c
# ╠═f15ebcd9-0755-4c9f-a417-21eb20e707c8
# ╠═cc5e59c1-f39f-40b6-b9b8-58d3adb8a54e
