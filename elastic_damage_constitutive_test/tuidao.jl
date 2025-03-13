### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 4eb112de-2fe9-4206-b08d-56c13a110731
using SymPy

# ╔═╡ 7ab25dea-8d8f-49be-a7af-9439a20cd37f
d, ν = Sym("d, ν")

# ╔═╡ 9875a55a-6bc9-4cae-bd1f-64f5b51263fb
function ×(A, B)
	return[
		A[1]*B[1] + 2A[6]*B[5]
		A[2]*B[2] + 2A[5]*B[6]
		A[3]*B[3]
		A[4]*B[4]
		A[5]*B[1] + A[2]*B[5]
		A[6]*B[2] + A[1]*B[6]
	]
end

# ╔═╡ e5ba3b71-044b-4f00-a47a-96e3062911c1
function ⦼(A)
	l = A[1] * A[2] - 2A[5] * A[6]
	return[
		A[2] / l
		A[1] / l
		1 / A[3]
		1 / A[4]
		-A[5] / l
		-A[6] / l
	]
end

# ╔═╡ 258b612d-f907-49d5-9346-7b4fd5648954
I = [1, 1, 1, 1, 0, 0]

# ╔═╡ 819d79a1-36e0-42f7-b50b-b6783a15ecbf
⦼(I)

# ╔═╡ 4ae325e3-fa0b-490d-a6cd-b6fab5dc2538
@syms a,b,c,e,f,cn, ω, E, μ12, μ13, k

# ╔═╡ 5d37c7c2-3d9a-4ad3-bf3e-35fa5ae3f444
En = k * E

# ╔═╡ d37d69d5-46fb-42bd-93e0-84841a660460
Shom = [(1-ν)/E, 1/En, (1+ν)/E, 1/(μ12+μ13), -ν/E, -ν/E]

# ╔═╡ 694520dc-d8bb-42ea-845e-b41a0089dca4
Chom = ⦼(Shom)

# ╔═╡ 5c24aa60-b7fb-46f0-9d24-fff07128b9a2
simplify.(Chom)

# ╔═╡ 5df673ea-5d21-4cf4-9bec-0c14753796f7
C0 = [1/ν, (1-ν)/ν, (1-2ν)/ν, (1-2ν)/ν, 1, 1] * E * ν / (1+ν) / (1-2ν)

# ╔═╡ f1a99390-e372-43c1-aa26-1f7795ceb0d6
r = 16(1-ν) / 3(1-2ν)

# ╔═╡ 91f9d652-0d8e-429f-ab34-1d323f52c4d8
ke = simplify((1 - r * (1 - ν) * d) / (1 - 2r * ν^2 * d))

# ╔═╡ 1a7cd126-b97f-4ac3-a270-22f43d70ee35
string(ke)

# ╔═╡ f3c7e7ab-89cf-4922-9d75-cae7fc041277
(3-6*0.2)/16(1-0.2)^2

# ╔═╡ a4f39aa0-0ac2-4864-8645-fb467b2f1271
z1 = 1 - 16(1 - ν^2)d / ( 9(1 - 2ν) + 16//3*(1 + ν)^2 * d )

# ╔═╡ ce054abf-df3c-49a0-a1b0-31e473551101
x1 = simplify(diff(z1, d))

# ╔═╡ 6ea8cd84-3718-49dc-af2f-4fae6b13739a
string(x1)

# ╔═╡ 21e0701b-3aaf-4bf8-ab3a-ec841d19966b
y1 = simplify(diff(x1, d))

# ╔═╡ 3a187170-1d02-4224-9486-37dcfd3051dc
string(y1)

# ╔═╡ 7c1bdc92-44b4-4930-8c1c-d8e4dfa69e11
ϑ = (5 - ν) / 3

# ╔═╡ b861e05f-ca2b-4cfc-8693-3269997be593
z2 = 1 - 480(1 - ν)ϑ * d / ( 225(2 - ν) + 64(4 - 5ν)ϑ * d )

# ╔═╡ 4f584ec5-4d0c-4fb0-a212-87967eb39f6e
x2 = simplify(diff(z2, d))

# ╔═╡ 498cd658-ff08-4392-9d77-27dc6bbec7a2
string(x2)

# ╔═╡ 462b46eb-aa40-401e-8f75-3df505ecdfe0
y2 = simplify(diff(x2, d))

# ╔═╡ 9ac4c13f-c912-4af4-a0bc-6f2059cc3f2c
factor(y2)

# ╔═╡ a1a978a3-5dc1-480d-a500-b215057e898a
string(y2)

# ╔═╡ 445f06a5-8cad-4c63-9900-7b77c902ec31
z3 = 1 - 480(1 - ν)d / ( 225(2 - ν) + 64(4 - 5ν)d )

# ╔═╡ d4e22b71-1731-46ba-84bb-439615e42614
factor(z3)

# ╔═╡ f84fcc99-397e-441c-9506-63e0e0203bb0
x3 = factor(diff(z3, d))

# ╔═╡ 58bb6760-84e5-4e9c-8725-1f2aed653fac
string(x3)

# ╔═╡ cc82c68d-cbac-4f2e-aff7-01ce45e96d29
y3 = factor(diff(x3, d))

# ╔═╡ c2c29bc0-ac1d-4d9c-bbad-9976528e8de9
string(y3)

# ╔═╡ 63918bdc-bc8e-41c5-8db4-34b9441faad1
β1 = 16(1 - ν^2) / 9(1 - 2ν)

# ╔═╡ f3cdf174-2574-4de0-9494-d093c650b572
32(1 - ν) / 15(2 - ν)

# ╔═╡ c9c54154-7641-4b57-ad0d-5742bb95fe2a
factor(1 - β1 * d)

# ╔═╡ 16c51893-de33-4342-9ef3-fc8eacd6a983
factor(1 / (1 + β1 * d))

# ╔═╡ 1ffd9e5c-edda-48fd-b7fc-5137964bd1c8
factor(diff(1 / (1 + β1 * d), d))

# ╔═╡ dad21832-b5ec-47d1-8545-9529b946ad12
simplify(diff(1 / (1 + β1 * d), d))

# ╔═╡ 9eab6452-2cfc-4b8a-bdbf-24dbe3b42ea2
rc, dc = Sym("rc, dc")

# ╔═╡ c1827768-070c-470c-bdfc-a860864ac949
ξ = d / dc

# ╔═╡ 47f7334b-1b4c-4924-9e7c-c94d04575759
string(factor(rc * 4ξ / (1 + ξ)^2))

# ╔═╡ 864d0012-9f39-4868-89d2-ed6a59d43c2a
string(factor(diff(rc * 4ξ / (1 + ξ)^2, d)))

# ╔═╡ 107f496e-2a52-43d5-aeac-85e644223f81
factor(diff(rc * 4ξ / (1 + ξ)^2, d))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
SymPy = "~2.2.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.0-rc3"
manifest_format = "2.0"
project_hash = "8a75300aca07d42c55ef9694fd4ed906d8fac190"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CommonEq]]
git-tree-sha1 = "6b0f0354b8eb954cdba708fb262ef00ee7274468"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.1"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "b19db3927f0db4151cb86d073689f2428e524576"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.10.2"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "9816a3826b0ebf49ab4926e2b18842ad8b5c8f04"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.96.4"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.SymPy]]
deps = ["CommonEq", "CommonSolve", "LinearAlgebra", "PyCall", "SpecialFunctions", "SymPyCore"]
git-tree-sha1 = "d35b297be048dfac05bcff29e55d6106808e3c5a"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "2.2.0"

[[deps.SymPyCore]]
deps = ["CommonEq", "CommonSolve", "Latexify", "LinearAlgebra", "Markdown", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "bef92ec4c31804bdc9c44cb00eaf0348eac383fb"
uuid = "458b697b-88f0-4a86-b56b-78b75cfb3531"
version = "0.2.5"

    [deps.SymPyCore.extensions]
    SymPyCoreTermInterfaceExt = "TermInterface"

    [deps.SymPyCore.weakdeps]
    TermInterface = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─4eb112de-2fe9-4206-b08d-56c13a110731
# ╠═7ab25dea-8d8f-49be-a7af-9439a20cd37f
# ╠═9875a55a-6bc9-4cae-bd1f-64f5b51263fb
# ╠═e5ba3b71-044b-4f00-a47a-96e3062911c1
# ╠═258b612d-f907-49d5-9346-7b4fd5648954
# ╠═819d79a1-36e0-42f7-b50b-b6783a15ecbf
# ╠═4ae325e3-fa0b-490d-a6cd-b6fab5dc2538
# ╠═5d37c7c2-3d9a-4ad3-bf3e-35fa5ae3f444
# ╠═d37d69d5-46fb-42bd-93e0-84841a660460
# ╠═694520dc-d8bb-42ea-845e-b41a0089dca4
# ╠═5c24aa60-b7fb-46f0-9d24-fff07128b9a2
# ╠═5df673ea-5d21-4cf4-9bec-0c14753796f7
# ╠═f1a99390-e372-43c1-aa26-1f7795ceb0d6
# ╠═91f9d652-0d8e-429f-ab34-1d323f52c4d8
# ╠═1a7cd126-b97f-4ac3-a270-22f43d70ee35
# ╠═f3c7e7ab-89cf-4922-9d75-cae7fc041277
# ╠═a4f39aa0-0ac2-4864-8645-fb467b2f1271
# ╠═ce054abf-df3c-49a0-a1b0-31e473551101
# ╠═6ea8cd84-3718-49dc-af2f-4fae6b13739a
# ╠═21e0701b-3aaf-4bf8-ab3a-ec841d19966b
# ╠═3a187170-1d02-4224-9486-37dcfd3051dc
# ╠═7c1bdc92-44b4-4930-8c1c-d8e4dfa69e11
# ╠═b861e05f-ca2b-4cfc-8693-3269997be593
# ╠═4f584ec5-4d0c-4fb0-a212-87967eb39f6e
# ╠═498cd658-ff08-4392-9d77-27dc6bbec7a2
# ╠═462b46eb-aa40-401e-8f75-3df505ecdfe0
# ╠═9ac4c13f-c912-4af4-a0bc-6f2059cc3f2c
# ╠═a1a978a3-5dc1-480d-a500-b215057e898a
# ╠═445f06a5-8cad-4c63-9900-7b77c902ec31
# ╠═d4e22b71-1731-46ba-84bb-439615e42614
# ╠═f84fcc99-397e-441c-9506-63e0e0203bb0
# ╠═58bb6760-84e5-4e9c-8725-1f2aed653fac
# ╠═cc82c68d-cbac-4f2e-aff7-01ce45e96d29
# ╠═c2c29bc0-ac1d-4d9c-bbad-9976528e8de9
# ╠═63918bdc-bc8e-41c5-8db4-34b9441faad1
# ╠═f3cdf174-2574-4de0-9494-d093c650b572
# ╠═c9c54154-7641-4b57-ad0d-5742bb95fe2a
# ╠═16c51893-de33-4342-9ef3-fc8eacd6a983
# ╠═1ffd9e5c-edda-48fd-b7fc-5137964bd1c8
# ╠═dad21832-b5ec-47d1-8545-9529b946ad12
# ╠═9eab6452-2cfc-4b8a-bdbf-24dbe3b42ea2
# ╠═c1827768-070c-470c-bdfc-a860864ac949
# ╠═47f7334b-1b4c-4924-9e7c-c94d04575759
# ╠═864d0012-9f39-4868-89d2-ed6a59d43c2a
# ╠═107f496e-2a52-43d5-aeac-85e644223f81
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
