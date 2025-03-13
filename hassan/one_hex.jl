### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 54a214aa-6454-11ef-2d38-1713c55850ac
begin
	import Pkg
	Pkg.activate(temp=true)
	using FEM, DelimitedFiles, Plots
end

# ╔═╡ 5aed4d17-3c7e-4f25-b6f7-888c505e4633
md"""
    4----------3
    |\         |\
    | \        | \
    |  \       |  \
    |   8----------7
    |   |      |   |
    1---+------2   |
     \  |       \  |
      \ |        \ |
       \|         \|
        5----------6
"""

# ╔═╡ 2ef60cfe-6592-437b-ae44-d4dd6ac6d109
begin
nodes = [0 0 0.0
         1 0 0
         1 0 1
         0 0 1
         0 -1 0
         1 -1 0
         1 -1 1
         0 -1 1]'
mesh1 = fem_mesh_gen(nodes)
end

# ╔═╡ 20dee44a-5280-4c65-b6c8-a82466791694
begin
    E = 7e4
	ν = 0.15
	rc = 9.36e-3
    dc = 7.0
	l = 3.0 # length scale 单位mm
    cl = elastic_damage{Dim3}(:MT, E, ν, Rd_gen(1, rc, dc))
end

# ╔═╡ 319fb310-0b49-4a21-91b7-7ac1e363495f
begin
e = [1 2 3 4 5 6 7 8]'
fem1 = elastoplastic(H1, e, nodes, 3, 3, 3, E = E, ν = ν, nstatev = 1)
end

# ╔═╡ aba593c2-efdd-4ae1-b07e-db249525f792
rcd = fem_recorder("t2", "output", fem1.quadrature_points)

# ╔═╡ a6ab0161-999f-47f0-aa58-5754ca6bee7b
begin
fixed_dofs = [1, 2, 3, 5, 6, 15, 18]
load_dofs = [9, 12, 21, 24]
max_load = 4e-2
nincr = 100
load_per_step = max_load / nincr
load_value = fill(load_per_step, length(load_dofs))
solve1 = simple_analysis(fixed_dofs, load_dofs, load_value, mesh1.dofs, tolerance = 1e-6)
end

# ╔═╡ 09fdfe12-1cf3-4529-ba64-43ad58cd1ef6
displacement_control!(rcd, solve1, fem1, mesh1, cl, nincr = nincr)

# ╔═╡ fc73426c-d146-4ef4-bade-4e1d3410801e
rcd, mesh1, fem1

# ╔═╡ 296388c2-387b-4481-b579-f376b942c8a2
record!(rcd)

# ╔═╡ 5f3fbc95-fb8f-42fa-a9c5-5b59fe29af5a
plot2vtu(rcd, mesh1, fem1.element_table)

# ╔═╡ 4a5bedb2-c2aa-4dbb-9700-7ac62b02ad7f
plot2vtu(rcd, mesh1, 1:nincr, fem1.element_table)

# ╔═╡ 2befb85b-5565-47e8-ba84-dd8863b8a4f4
rcd

# ╔═╡ Cell order:
# ╠═54a214aa-6454-11ef-2d38-1713c55850ac
# ╟─5aed4d17-3c7e-4f25-b6f7-888c505e4633
# ╠═2ef60cfe-6592-437b-ae44-d4dd6ac6d109
# ╠═20dee44a-5280-4c65-b6c8-a82466791694
# ╠═319fb310-0b49-4a21-91b7-7ac1e363495f
# ╠═aba593c2-efdd-4ae1-b07e-db249525f792
# ╠═a6ab0161-999f-47f0-aa58-5754ca6bee7b
# ╠═09fdfe12-1cf3-4529-ba64-43ad58cd1ef6
# ╠═fc73426c-d146-4ef4-bade-4e1d3410801e
# ╠═296388c2-387b-4481-b579-f376b942c8a2
# ╠═5f3fbc95-fb8f-42fa-a9c5-5b59fe29af5a
# ╠═4a5bedb2-c2aa-4dbb-9700-7ac62b02ad7f
# ╠═2befb85b-5565-47e8-ba84-dd8863b8a4f4
