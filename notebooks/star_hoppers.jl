### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ 03688332-bcc4-4c84-8924-7f1610db077e
# imports
begin
	using Dierckx
	using Markdown
	using Measurements
	using Unitful
	using UnitfulAngles
	using UnitfulAstro
end

# ╔═╡ b1450a24-3120-11ec-2c56-d7b2f96fdf46
md"""
# Star Hoppers

This notebook will follow along with some derivations of [Kratter & Perets (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...753...91K).
"""

# ╔═╡ 2dbe3b2b-f33b-46f3-98f3-869c539451db
M_A = (2.063 ± 0.023)u"Msun"

# ╔═╡ 6eef835d-49b3-4ffe-bc3a-dc7d5938fe5b
M_B_WD = (1.018 ± 0.011)u"Msun"

# ╔═╡ c79b9551-4ccc-4ea4-8e41-1c56cd5e4c72
md"white dwarf initial-final mass relation (IFMR) from [Cummings et al. (eqn. 6; 2018)](https://ui.adsabs.harvard.edu/abs/2018ApJ...866...21C)"

# ╔═╡ 6062c398-de6a-457f-a7c5-e8ce13d7bd3e
function MIST_highmass_IFMR(M_WD)
	b = (0.471 ± 0.077)u"Msun"
	m = (0.107 ± 0.016)
	return (M_WD - b) / m
end

# ╔═╡ eb6f04d6-596f-497e-93ab-fdcde0e92d8d
M_B_MS = MIST_highmass_IFMR(M_B_WD)

# ╔═╡ df958bae-32e7-47b6-8e9f-95698151113c
md"""
Adiabatic expansion factor

```math
f = \frac{a_f}{a_i} = \frac{M_\mathrm{MS}}{M_\mathrm{WD}}
```
"""

# ╔═╡ f6d29ca4-b536-49a9-9391-81b44d3e7413
f = M_B_MS / M_B_WD

# ╔═╡ c431232c-b181-483a-95a7-25426d865281


# ╔═╡ 8b14cc60-5373-4ca9-bc90-932e1a729e0e
# semimajor axis from GAIA parallax and Bond et al. 2017
begin
	plx = (0.37449 ± 0.00023)
	a_arc = (7.4957 ± 0.0025)
	a = (a_arc / plx) * u"AU"
end

# ╔═╡ b5d25245-2c5c-41ce-a255-8fb5aa7b4ee8
(M_B_WD + M_A) / (M_B_MS + M_A) |> inv

# ╔═╡ 4c3173e9-3599-4bba-864d-354c369a2a3e
a_initial = a * (M_B_WD + M_A) / (M_B_MS + M_A)

# ╔═╡ 83e74886-9204-49ae-8f65-e9550e442c39
@. a_initial * (1 + [0.59, -0.59])

# ╔═╡ a5527e3e-15f1-4f14-9c2a-f77706dfda3e
md"""
Forbidden region defined by minimum separation required to (1) escape engulfment and (2) escape tidal shredding. The first is defined by the max giant branch (GB) radius, and the second is dependent on the companion mass, but for giant planets we assume some flat percentage extra (here, 20% is based off Veras 2016)
"""

# ╔═╡ 2f6a3dfc-ca25-455e-ad61-df8eaa049e08
max_GB_radius = (5.104 ± 0.075)u"AU"

# ╔═╡ 677ad47b-538a-41e4-844b-8fff1a3d6314
tide_factor = 0.2

# ╔═╡ 537a9087-ddb1-45e8-a872-593f6fad93da
stable_radius(M1, M2, a) = (0.464 - 0.38 * M2 / (M1 + M2)) * a	

# ╔═╡ 11d5f746-56c5-49dd-971a-fd90e47a7952
R_P = stable_radius(M_B_MS, M_A, a_initial)

# ╔═╡ 1ac684d7-f0c6-426f-8318-075352c92ca0
R_S = stable_radius(M_A, M_B_MS, a_initial)

# ╔═╡ 88f576cf-696d-4d7f-8123-729bdef2d9a6
area_ratio = R_S^2 / (R_P^2 + R_S^2)

# ╔═╡ 2c15891e-c3bf-4a43-b4b1-4567ad7c5ed9
coll_prob_m1 = Spline1D(
	[0.5, 0.55, 0.6, 0.65, 0.7],
	[0.4, 0.49, 0.63, 0.75, 0.8],
	k = 1,
	bc="extrapolate"
)

# ╔═╡ 28f682d4-f9bd-4a00-a2ec-c6be9ceaeb46
coll_prob_m2 = Spline1D(
	[0.5, 0.55, 0.6, 0.65, 0.7],
	[0.53, 0.5, 0.27, 0.22, 0.17],
	k = 1,
	bc="extrapolate"
)

# ╔═╡ a365cdca-14fb-4889-afde-431f8b292ba2
@uncertain coll_prob_m1(area_ratio)

# ╔═╡ 2a550342-7b8c-4dde-a89f-e7e7e63d1bc1
@uncertain coll_prob_m2(area_ratio)

# ╔═╡ 29b9a599-a519-4b72-afd1-3aa1d872a8cf
min_escape_separation_pre = max_GB_radius * (1 + tide_factor)

# ╔═╡ c61024e3-769b-4968-a1ba-a5509de1848a
min_escape_separation_post = min_escape_separation_pre * f

# ╔═╡ c6de59a9-57b2-40ff-8732-d4a01f11db17
(M_B_MS + M_A) / 2.5u"Msun" * ustrip(u"Msun", M_A)^(3/2) / sqrt(1.711) / (a_initial / 100u"AU")^(3/2)

# ╔═╡ 806d35ba-1efb-44eb-9153-64f6992bfe31
(M_B_WD + M_A) / 2.5u"Msun" * ustrip(u"Msun", M_A)^(3/2) / sqrt(1.711) / (a / 100u"AU")^(3/2)

# ╔═╡ 286abf89-2ffc-42e8-8f89-b8b58ffe507e
begin
	μ = M_B_MS / (M_A + M_B_MS)
	eb = 0.59
	ac_ab = (0.464 ± 0.006) + (-0.38 ± 0.01) * μ +
		(-0.631 ± 0.034) * eb + (0.586 ± 0.061) * eb * μ +
		(0.15 ± 0.041) * eb^2 + (-0.198 ± 0.047) * μ * eb^2
	ac = ac_ab * a
end

# ╔═╡ Cell order:
# ╠═03688332-bcc4-4c84-8924-7f1610db077e
# ╟─b1450a24-3120-11ec-2c56-d7b2f96fdf46
# ╟─2dbe3b2b-f33b-46f3-98f3-869c539451db
# ╟─6eef835d-49b3-4ffe-bc3a-dc7d5938fe5b
# ╟─c79b9551-4ccc-4ea4-8e41-1c56cd5e4c72
# ╠═6062c398-de6a-457f-a7c5-e8ce13d7bd3e
# ╠═eb6f04d6-596f-497e-93ab-fdcde0e92d8d
# ╟─df958bae-32e7-47b6-8e9f-95698151113c
# ╠═f6d29ca4-b536-49a9-9391-81b44d3e7413
# ╠═c431232c-b181-483a-95a7-25426d865281
# ╠═8b14cc60-5373-4ca9-bc90-932e1a729e0e
# ╠═b5d25245-2c5c-41ce-a255-8fb5aa7b4ee8
# ╠═4c3173e9-3599-4bba-864d-354c369a2a3e
# ╠═83e74886-9204-49ae-8f65-e9550e442c39
# ╟─a5527e3e-15f1-4f14-9c2a-f77706dfda3e
# ╟─2f6a3dfc-ca25-455e-ad61-df8eaa049e08
# ╟─677ad47b-538a-41e4-844b-8fff1a3d6314
# ╠═537a9087-ddb1-45e8-a872-593f6fad93da
# ╠═11d5f746-56c5-49dd-971a-fd90e47a7952
# ╠═1ac684d7-f0c6-426f-8318-075352c92ca0
# ╠═88f576cf-696d-4d7f-8123-729bdef2d9a6
# ╠═2c15891e-c3bf-4a43-b4b1-4567ad7c5ed9
# ╠═28f682d4-f9bd-4a00-a2ec-c6be9ceaeb46
# ╠═a365cdca-14fb-4889-afde-431f8b292ba2
# ╠═2a550342-7b8c-4dde-a89f-e7e7e63d1bc1
# ╠═29b9a599-a519-4b72-afd1-3aa1d872a8cf
# ╠═c61024e3-769b-4968-a1ba-a5509de1848a
# ╠═c6de59a9-57b2-40ff-8732-d4a01f11db17
# ╠═806d35ba-1efb-44eb-9153-64f6992bfe31
# ╠═286abf89-2ffc-42e8-8f89-b8b58ffe507e
