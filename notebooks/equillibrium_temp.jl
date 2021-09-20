### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 7927ff00-1806-11ec-2787-4db3c2bd8448
# imports
begin
	using Unitful, UnitfulAstro
	using UnitfulRecipes
	using Revise
	using LinearAlgebra
	using Pkg
	Pkg.develop(path="/Users/miles/dev/julia/Transits")
	using Transits
	using Plots
	using Markdown
end

# ╔═╡ 87591ebf-f36b-4643-bbf1-a961ed1825ff
md"""
# Equillibrium temperature in Sirius

When considering a planet orbiting Sirius B, it is interesting to consider the irradiation from Sirius A. This is especially significant due to Sirius B's eccentric orbit (``e=0.6``), leading to an overall range of separations from 7.9-32 AU.

The calculation for equilibrium temperature from a single star is
```math
T_\mathrm{eq} = \left[\frac{L(1 - A_B)}{16\sigma\pi a^2} \right]^{1/4}
```

we assume that adding an additional star will linearly add to the incident flux, so we can extend this to two stars with
```math
T_\mathrm{eq} = \sum_i{\left[\frac{L_i(1 - A_B)}{16\sigma\pi a_i^2} \right]^{1/4}}
```

This does not take into effect any reheating or geometric phase effects, which means this is likely an *upper limit*.
"""

# ╔═╡ b4375a6e-c669-4f64-aa9c-72375bf96868
begin
	plx = 0.3789 # arcsec
	dist = inv(plx)u"pc" # parsec
end

# ╔═╡ 0bdbeded-0d6f-426b-80bd-0d20d058baf8
# don't include i, Ω, or ω because we don't care about the 
# line-of-sight rotation, only the relative position
# orbital elements from Bond et al. 2017
orbit = KeplerianOrbit(
	period=50.1284u"yr",
	a=7.4957u"AU" / plx, # convert arcsec to AU
	#incl=136.336u"°",
	#Omega=45.400u"°",
	t_0=1994.5715u"yr",
	ecc=0.59142,
	omega=0u"°"
	#omega=149.161u"°"
)

# ╔═╡ 0dd1df5c-7faa-4795-b4cd-0ebb9ba9e0c0
function relative_position(orbit, t)
	# orbit is in xzy because it wasn't rotated
	xzy = Transits.Orbits.relative_position(orbit, t)
	return [xzy[1], xzy[3]]
end

# ╔═╡ 181ea5bd-18f0-47ae-8911-9c823b8b5531
# plot orbit
begin
	t = range(0u"d", orbit.period, length=1000)
	t_today = 2020u"yr" + 260u"d" |> u"d"
	# calculate separations
	pos = mapreduce(ti -> relative_position(orbit, ti), hcat, t)
	plot(eachrow(pos)..., 
		c=:black, 
		leg=false, 
		xlabel="Δx [AU]", ylabel="Δy [AU]")
	scatter!([0], [0], m=:star, ms=10, c=:black)
	pos_today = relative_position(orbit, t_today)
	scatter!(pos_today[1:1], pos_today[2:2], m=:star, ms=5, c=:black)
end

# ╔═╡ 1235c53f-e01f-457d-9057-7587af56a9f8
function eq_temp(orbit, t, dx, dy)
	pos = relative_position(orbit, t)u"AU"
	d = [dx, dy]
	a1 = norm(d)        # dist to Sirius A
	a2 = norm(d .- pos) # dist to Sirius B
	# assume albedo is negligible
	L1 = 25.4u"Lsun"
	L2 = 0.056u"Lsun"
	Ts = (L1 / (16π * Unitful.σ * a1^2))^(1/4) + 
	     (L2 / (16π * Unitful.σ * a2^2))^(1/4)
	return Ts |> u"K"
end

# ╔═╡ 7738d09a-b7c5-4a38-b0fe-aafcba7e092a
begin
	dxs = range(-10, 35, length=500)u"AU"
	dys = range(-17, 17, length=500)u"AU"
end

# ╔═╡ 0dfa9c8b-0dd1-4967-99ce-6e729a3506a9
function eq_temp(orbit, t, dist)
	pos = relative_position(orbit, t)u"AU"
	d = norm(pos) # total orbital distance
	a1 = d - dist # distance to Sirius A
	a2 = dist     # distance to Sirius B
	# assume albedo is negligible
	L1 = 25.4u"Lsun"
	L2 = 0.056u"Lsun"
	T1 = (L1 / (16π * Unitful.σ * a1^2))^(1/4) 
	T2 = (L2 / (16π * Unitful.σ * a2^2))^(1/4)
	return (T1 + T2, T1, T2) .|> u"K"
end

# ╔═╡ 57a6fbf8-0b5b-4a0e-b60a-e7d367692446
let temps
	temps = eq_temp.(orbit, t_today, dxs', dys)
	heatmap(dxs, dys, log10.(ustrip.(temps)), 
		# levels=100,
		xlim=extrema(dxs),
		ylim=extrema(dys),
		c=:acton,
		cbar_title="log10(Teq)",
		title="today"
	)
	plot!(eachrow(pos)..., c=:white, leg=false, xlabel="Δx [AU]", ylabel="Δy [AU]")
	scatter!([0], [0], m=:star, ms=10, c=:white)
	scatter!(pos_today[1:1], m=:star, ms=5, pos_today[2:2], c=:white)
end

# ╔═╡ 56cd4737-3bdc-43da-a1fc-530baca135a3
let temps
	t_per = orbit.t_0 + orbit.t_ref
	temps = eq_temp.(orbit, t_per, dxs', dys)
	heatmap(dxs, dys, log10.(ustrip.(temps)), 
		# levels=100,
		xlim=extrema(dxs),
		ylim=extrema(dys),
		c=:acton,
		cbar_title="log10(Teq)",
		title="periapsis"
	)
	plot!(eachrow(pos)..., c=:white, leg=false, xlabel="Δx [AU]", ylabel="Δy [AU]")
	scatter!([0], [0], m=:star, ms=10, c=:white)
	pos_per = relative_position(orbit, t_per)
	scatter!(pos_per[1:1], pos_per[2:2], m=:star, ms=5, c=:white)
end

# ╔═╡ 0523a4b3-7ecc-444e-853f-3adb46b8c04d
let temps
	t_apo = orbit.t_0 + orbit.t_ref + orbit.period/2
	temps = eq_temp.(orbit, t_apo, dxs', dys)
	heatmap(dxs, dys, log10.(ustrip.(temps)), 
		# levels=100,
		xlim=extrema(dxs),
		ylim=extrema(dys),
		c=:acton,
		cbar_title="log10(Teq)",
		title="apoapsis"
	)
	plot!(eachrow(pos)..., c=:white, leg=false, xlabel="Δx [AU]", ylabel="Δy [AU]")
	scatter!([0], [0], m=:star, ms=10, c=:white)
	pos_apo = relative_position(orbit, t_apo)
	scatter!(pos_apo[1:1], pos_apo[2:2], m=:star, ms=5, c=:white)
end

# ╔═╡ cb972d70-1652-464a-9b35-215921674735
begin
	dists = range(0, 2, length=100)u"AU"
	plot(dists, map(d -> eq_temp(orbit, t_today, d)[1], dists), 
		lab="",
		ylabel="Teq",
		xlabel="separation from Sirius B",
		title="today",
		xlim=extrema(dists)
	)
	plot!(dists, map(d -> eq_temp(orbit, t_today, d)[2], dists),
		ls=:dash,
		lab="Sirius A"
	)
	plot!(dists, map(d -> eq_temp(orbit, t_today, d)[3], dists),
		ls=:dash,
		lab="Sirius B"
	)
end

# ╔═╡ 6c63961a-f6b6-47e7-8dad-854ff1746aeb
begin
	t_per = orbit.t_0 + orbit.t_ref
	plot(dists, map(d -> eq_temp(orbit, t_per, d)[1], dists), 
		lab="",
		ylabel="Teq",
		xlabel="separation from Sirius B",
		title="periapsis",
		xlim=extrema(dists)
	)
	plot!(dists, map(d -> eq_temp(orbit, t_per, d)[2], dists),
		ls=:dash,
		lab="Sirius A"
	)
	plot!(dists, map(d -> eq_temp(orbit, t_per, d)[3], dists),
		ls=:dash,
		lab="Sirius B"
	)
end

# ╔═╡ 860fc3e3-08d5-4cd5-8a4a-3b9db25a2b0e
begin
	t_apo = orbit.t_0 + orbit.t_ref + orbit.period/2
	plot(dists, map(d -> eq_temp(orbit, t_apo, d)[1], dists),
		lab="",
		ylabel="Teq",
		xlabel="separation from Sirius B",
		title="apoapsis",
		xlim=extrema(dists)
	)
	plot!(dists, map(d -> eq_temp(orbit, t_apo, d)[2], dists),
		ls=:dash,
		lab="Sirius A"
	)
	plot!(dists, map(d -> eq_temp(orbit, t_apo, d)[3], dists),
		ls=:dash,
		lab="Sirius B"
	)
end

# ╔═╡ d97e3642-f701-46e2-8c59-13148d2bebae
begin
	curves = hcat(
		map(d -> eq_temp(orbit, t_today, d)[1], dists),
		map(d -> eq_temp(orbit, t_per, d)[1], dists),
		map(d -> eq_temp(orbit, t_apo, d)[1], dists),
	)
	plot(dists, curves, 
		lab=["today" "periapsis" "apoapsis"], 
		ylabel="Teq", 
		xlabel="separation from Sirius B",
		xlim=extrema(dists)
	)
end

# ╔═╡ Cell order:
# ╟─87591ebf-f36b-4643-bbf1-a961ed1825ff
# ╠═b4375a6e-c669-4f64-aa9c-72375bf96868
# ╠═0bdbeded-0d6f-426b-80bd-0d20d058baf8
# ╠═0dd1df5c-7faa-4795-b4cd-0ebb9ba9e0c0
# ╟─181ea5bd-18f0-47ae-8911-9c823b8b5531
# ╠═1235c53f-e01f-457d-9057-7587af56a9f8
# ╠═7738d09a-b7c5-4a38-b0fe-aafcba7e092a
# ╟─57a6fbf8-0b5b-4a0e-b60a-e7d367692446
# ╟─56cd4737-3bdc-43da-a1fc-530baca135a3
# ╟─0523a4b3-7ecc-444e-853f-3adb46b8c04d
# ╠═0dfa9c8b-0dd1-4967-99ce-6e729a3506a9
# ╟─cb972d70-1652-464a-9b35-215921674735
# ╟─6c63961a-f6b6-47e7-8dad-854ff1746aeb
# ╟─860fc3e3-08d5-4cd5-8a4a-3b9db25a2b0e
# ╟─d97e3642-f701-46e2-8c59-13148d2bebae
# ╠═7927ff00-1806-11ec-2787-4db3c2bd8448
