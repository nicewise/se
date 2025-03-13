### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ daf13286-5fc0-11ef-16a4-2f0b4249d254
begin
	import Pkg
	Pkg.activate(temp=true)
	using FEM, DelimitedFiles, Plots
end

# ╔═╡ a3d61ed6-b5dc-431b-b721-4ed18d1e01d6
nincr = 100

# ╔═╡ 53131c09-e63c-40c7-94fb-5e15ef80652b
# ╠═╡ disabled = true
#=╠═╡
gradient_displacement_control!(rc, solve1, fem1, mesh1, cl, nincr = nincr)
  ╠═╡ =#

# ╔═╡ 4308759b-b779-4de5-8087-6d78d43b0f7c
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
	#cl = gradient_elastic_damage{Dim3}(:MT, E, ν, l, Rd_gen(1, rc, dc))
    cl = elastic_damage{Dim3}(:MT, E, ν, Rd_gen(1, rc, dc))

    mesh1 = fem_mesh_gen(nodes)
    fem1 = elastoplastic(H1, e, nodes, 3, 3, 3, E = E, ν = ν, nstatev = 2)
    rcd = fem_recorder("t2", "output_gradient", fem1.quadrature_points)

    top = findall(nodes[2, :].==maximum(nodes[2, :]))
    bottom = findall(nodes[2, :].==minimum(nodes[2, :]))
    fixed_dofs = vcat(3bottom.-2, 3bottom.-1, 3bottom) # 底部固定
    load_dofs = 3top.-1 # 顶部y方向加载
    max_load = 4e-2
    load_per_step = max_load / nincr
    load_value = fill(load_per_step, length(load_dofs))

    solve1 = simple_analysis(fixed_dofs, load_dofs, load_value, mesh1.dofs, tolerance = 1e-6)
    return rcd, cl, mesh1, fem1, solve1
end

# ╔═╡ 9c6c08bd-96a2-4004-ac74-3050f098cb7e
rcd, cl, mesh1, fem1, solve1 = model_gen(nincr)

# ╔═╡ b61793eb-5f1d-4c28-af99-96dbd2744d40
scatter(mesh1.nodes[1, :], mesh1.nodes[2, :], legend=false, axis=false)

# ╔═╡ 4781a19a-2124-4aeb-b809-fd9f8d972f75
displacement_control!(rcd, solve1, fem1, mesh1, cl, nincr = nincr)

# ╔═╡ 008678b4-2a0a-4eac-a8ea-a8c8cec1ac29
plot2vtu(rcd, mesh1, 1:nincr, fem1.element_table)

# ╔═╡ Cell order:
# ╟─daf13286-5fc0-11ef-16a4-2f0b4249d254
# ╠═b61793eb-5f1d-4c28-af99-96dbd2744d40
# ╠═a3d61ed6-b5dc-431b-b721-4ed18d1e01d6
# ╠═9c6c08bd-96a2-4004-ac74-3050f098cb7e
# ╠═53131c09-e63c-40c7-94fb-5e15ef80652b
# ╠═4781a19a-2124-4aeb-b809-fd9f8d972f75
# ╠═4308759b-b779-4de5-8087-6d78d43b0f7c
# ╠═008678b4-2a0a-4eac-a8ea-a8c8cec1ac29
