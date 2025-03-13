using FEM, DelimitedFiles
using Debugger

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
    fem1 = elastoplastic(H1, e, nodes, 3, 3, 3, E = E, ν = ν, nstatev = 1)
    rcd = fem_recorder("hassan1", "output", fem1.quadrature_points)

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

function analysis(nincr)
    rcd, cl, mesh1, fem1, solve1 = model_gen(nincr)
    displacement_control!(rcd, solve1, fem1, mesh1, cl, nincr = nincr)
    record!(rcd)
    return rcd, mesh1, fem1
end
function output(nincr, rcd, mesh, fe)
    for i in 1:nincr
        plot2vtu(rcd, mesh, i, fe.element_table)
        println(i, " done.")
    end
end
