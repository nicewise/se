### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 3de542fb-b092-493b-b236-0c5df62ab352
begin
	import Pkg
	Pkg.activate(temp=true)
end

# ╔═╡ d649e9f9-be03-4b76-90e1-0879639e8437
using FEM, DelimitedFiles, Plots

# ╔═╡ cdd4c5ad-1fd0-4379-b273-b1b76c0edd9d
abstract type bb <: AbstractBond end # bb .EQ. bond based

# ╔═╡ a54733d9-53c6-445d-9b2e-dffc1ffcb025
abstract type xbb <: AbstractBond end # xbb .EQ. extended bond based

# ╔═╡ 5b5d1630-5848-465e-ab5a-310b5354d6f3
abstract type AbstractPD{C<:AbstractBond, D<:AbstractDimension} end

# ╔═╡ d52127a2-c258-4159-b44d-7cbc85196098
struct bond_based{C<:AbstractBond, D<:AbstractDimension} <: AbstractPD{C, D}
    bond_table::Matrix{Int}
    node_index::Vector{Int}
    node_map::Dict{Int, Int}
    pdnodes::Vector{AbstractPoint{C, D}}
    bonds::Vector{AbstractCell{C, D}}
end

# ╔═╡ 80883113-a863-4471-a47f-08878437fc78
struct extended_bond_based{C<:AbstractBond, D<:AbstractDimension} <: AbstractPD{C, D}
    bond_table::Matrix{Int}
    node_index::Vector{Int}
    region::Vector{Vector{Int}}
    node_map::Dict{Int, Int}
    pdnodes::Vector{AbstractPoint{C, D}}
    bonds::Vector{AbstractCell{C, D}}
end

# ╔═╡ 9bcb047d-c47f-4ebe-ab99-343fa2269bb5


# ╔═╡ 4141cfab-0a67-4a07-ac38-b19cb3591345


# ╔═╡ 9fce762d-f608-4042-8580-0cc53a44d8e1


# ╔═╡ c63b0f08-b872-4ce6-b206-10a50d96af87


# ╔═╡ d48b8f2a-732f-48ba-a228-25ec71207b96


# ╔═╡ 76906d74-5d4e-40cb-a235-8487b0d98f0b


# ╔═╡ 45b0d159-3674-4634-8144-8d30155bf165


# ╔═╡ 4f723028-3fa4-481d-bd89-9187b8fb1909


# ╔═╡ 06be03ec-4a6b-4d96-875c-2e772a32fd9e
struct pd_mesh{D<:AbstractDimension} <: AbstractMesh{D}
    nodes::Matrix{Float64}
    dofs::Int
    inter_crack::Function
    thick::Float64
    δ::Float64
    dx::Float64
end

# ╔═╡ 5ea83478-f3a9-444c-b856-5abdd4f133fd
function create_inter_crack(crack::Matrix{Float64})
    if isempty(crack)
        return (n1::Vector{Float64}, n2::Vector{Float64}) -> true
    end

    return (n1::Vector{Float64}, n2::Vector{Float64}) -> begin
        ax, ay = n1[1], n1[2]
        bx, by = n2[1], n2[2]

        for k in 1:size(crack, 1)
            cx, cy = crack[k, 1], crack[k, 2]
            dx, dy = crack[k, 3], crack[k, 4]

            u = (cx - ax) * (by - ay) - (bx - ax) * (cy - ay)
            v = (dx - ax) * (by - ay) - (bx - ax) * (dy - ay)
            w = (ax - cx) * (dy - cy) - (dx - cx) * (ay - cy)
            z = (bx - cx) * (dy - cy) - (dx - cx) * (by - cy)

            if u * v < 0 && w * z < 0
                return false
            end
        end

        return true
    end
end

# ╔═╡ 13711a31-8289-48e5-97c3-aef269116c64
function pd_mesh_gen(nodes; thick = nothing, delta = nothing, dx = nothing, crack = zeros(0, 4))
    dimension, node_number = size(nodes) # 2xnnodes或3xnnodes
    # 维度检查和厚度验证
    @assert !(dimension == 2 && isnothing(thick)) "For 2D problem, `thick` must be provided."
    thick = dimension == 2 ? thick : 1.0
	@assert !isnothing(dx) "uniform mesh size `dx` must be provided."
	if isnothing(delta)
		delta = 1.5 * dx
	end
    dofs = dimension * node_number
	D = dimension == 1 ? Dim1 :
        dimension == 2 ? Dim2 :
        dimension == 3 ? Dim3 :
        nothing
	inter_crack = create_inter_crack(crack)
    pd_mesh{D}(nodes, dofs, inter_crack, thick, delta, dx)
end

# ╔═╡ 61ee7bad-f651-4c56-8dc0-ba137c53bfdc
elements = readdlm("elements", Int)'

# ╔═╡ d0787781-6db4-486f-843a-47bbe75f6732
nodes = readdlm("nodes")'

# ╔═╡ d625751d-39e9-4441-bb1d-f8231a4a5096
t = pd_mesh_gen(nodes, thick=1., dx=1)

# ╔═╡ af64f219-acf7-48e6-8868-367ad3ffdbd6
abstract type AbstractFailJudge end

# ╔═╡ a086e6a0-11e1-49c2-9507-fe24a13fa60f
struct critical_strecth <: AbstractFailJudge
	sc::Float64 # critical stretch
	max_fail_number::Int
end

# ╔═╡ ec1c28f4-9472-49a3-8b8c-8d97d47f650d
function bond_fail_judge!(pd::pd, fj::critical_strecth) # fj .EQ. fail judge
	# 创建包含键序号和伸长率的矩阵
    bonds = [i for i in 1:pd.bond_number]
    # 找到未断键的键
    valid_bonds = bonds[pd.bond_status]
    valid_stretch = pd.bond_stretch[pd.bond_status]

    # 找到可能断键的键
    potential_failures = valid_bonds[valid_stretch .> fj.sc]

    # 根据可能断键的数量确定断键序号
    if length(potential_failures) > fj.max_fail_number
        # 按伸长率降序排序，取前 max_fail_number 个键
	    fail_index = potential_failures[sortperm(valid_stretch[valid_stretch .> fj.sc], rev=true)[1:fj.max_fail_number]]
    else
	    fail_index = potential_failures
    end

    # 更新断键键的连接状态为 false
    pd.bond_status[fail_index] .= false
    return fail_index
end

# ╔═╡ 341ee260-29d6-4dae-9184-2cf393689583
function build_fem_model()
    elements = readdlm("elements", Int)'
    nodes = readdlm("nodes")' # 单位mm

    E = 3e4 # 单位 MPa
    ν = 1/3
    thick = 1 # 单位 mm

    mesh1 = fem_mesh_gen(nodes, thick = thick)
    fem1 = linear_elastic(q1, elements, 2, 2, E = E, ν = ν, isplanestrain=false)

    top = findall(nodes[2, :].==maximum(nodes[2, :]))
    bottom = findall(nodes[2, :].==minimum(nodes[2, :]))
    fixed_dofs = vcat(2top, 2bottom.-1, 2bottom)
    load_dofs = 2top.-1
    max_load = 2e-2
    load_per_step = 2e-2
    max_load_factor = max_load / load_per_step
    load_value = fill(load_per_step, length(load_dofs))

    solve1 = simple_analysis(fixed_dofs, load_dofs, load_value, mesh1.dofs, tolerance = 1e-6)
    return mesh1, fem1, solve1
end

# ╔═╡ 936171d9-f325-46d0-a4e8-de041ae96c15
mesh1, fem1, solve1 = build_fem_model()

# ╔═╡ Cell order:
# ╠═3de542fb-b092-493b-b236-0c5df62ab352
# ╠═d649e9f9-be03-4b76-90e1-0879639e8437
# ╠═cdd4c5ad-1fd0-4379-b273-b1b76c0edd9d
# ╠═a54733d9-53c6-445d-9b2e-dffc1ffcb025
# ╠═5b5d1630-5848-465e-ab5a-310b5354d6f3
# ╠═d52127a2-c258-4159-b44d-7cbc85196098
# ╠═80883113-a863-4471-a47f-08878437fc78
# ╠═9bcb047d-c47f-4ebe-ab99-343fa2269bb5
# ╠═4141cfab-0a67-4a07-ac38-b19cb3591345
# ╠═9fce762d-f608-4042-8580-0cc53a44d8e1
# ╠═c63b0f08-b872-4ce6-b206-10a50d96af87
# ╠═d48b8f2a-732f-48ba-a228-25ec71207b96
# ╠═76906d74-5d4e-40cb-a235-8487b0d98f0b
# ╠═45b0d159-3674-4634-8144-8d30155bf165
# ╠═4f723028-3fa4-481d-bd89-9187b8fb1909
# ╠═06be03ec-4a6b-4d96-875c-2e772a32fd9e
# ╠═13711a31-8289-48e5-97c3-aef269116c64
# ╠═5ea83478-f3a9-444c-b856-5abdd4f133fd
# ╠═61ee7bad-f651-4c56-8dc0-ba137c53bfdc
# ╠═d0787781-6db4-486f-843a-47bbe75f6732
# ╠═d625751d-39e9-4441-bb1d-f8231a4a5096
# ╠═af64f219-acf7-48e6-8868-367ad3ffdbd6
# ╠═a086e6a0-11e1-49c2-9507-fe24a13fa60f
# ╠═ec1c28f4-9472-49a3-8b8c-8d97d47f650d
# ╠═341ee260-29d6-4dae-9184-2cf393689583
# ╠═936171d9-f325-46d0-a4e8-de041ae96c15
