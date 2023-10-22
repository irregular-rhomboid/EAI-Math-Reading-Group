### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ b51af2cc-653e-11ee-19cc-bb1f00f82e59
using Pkg; Pkg.activate()

# ╔═╡ 6d86c249-016e-4ac0-ae42-5ee5afa305d9
begin
	using LinearAlgebra
	using Plots
	#using PlotlyJS # You'll need to install it to run the notebook but after that Plots will handle things
	using StatsPlots
	using Distributions
	using Statistics
	using PlutoUI
	using LaTeXStrings
	using Random
end

# ╔═╡ 86d7fbae-2b06-4506-b635-0f3ce5cc4c4e
md"# Random Matrix Theory: Wishart Matrices and the Marcenko-Pastur distribution"

# ╔═╡ 47183ef7-a6a1-4626-8449-12950733c6cc
plotlyjs()

# ╔═╡ 84849894-e7a0-4698-a8e7-9ab08db88838
begin
	xmin = -2.0
	xmax = 4.0
	ymin = -2.0
	ymax = 2.0
	step = 0.01
	xs = xmin:step:xmax
	ys = ymin:step:ymax
end

# ╔═╡ de3d0390-844e-4053-996c-c83349b82503
λ₋(q) = (1-sqrt(q))^2

# ╔═╡ c2b60494-6de3-4e31-8dbb-911457d818e6
λ₊(q) = (1+sqrt(q))^2

# ╔═╡ e8aad398-2dff-422f-b754-2e35dea55019
g(z,q) = (z - (1-q) - sqrt((z-(1-sqrt(q))^2))*sqrt((z-(1+sqrt(q))^2))) / (2*q*z)

# ╔═╡ 67013fbf-8b8b-4c8c-a395-9152aa760054
@bind q Slider(0.0:0.01:2.0, default=1.0, show_value=true)

# ╔═╡ 62e2f4ab-9d93-49b3-a0b8-8c28156beb86
g_(x,y) = g(x+y*im, q)

# ╔═╡ 90ff0901-8b0e-4306-b9c5-5d5bb7b2872a
begin
	#plotlyjs()
	local z = (x,y) -> min(4, abs(g_(x,y)))
	local p = Plots.surface(xs,ys, z,
		title="Stieltjes Transform of a White Wishart Matrix",
		size=(700,600),
		xlabel="Re(z)",
		ylabel="Im(z)",
		zlabel="|g(z)|"
	)
	p
end

# ╔═╡ 2ce15f2f-6e14-4a80-be68-acc92c7ea43f
begin
	#plotlyjs()
	local z = (x,y) -> min(10, abs(g_(x,y)))
	local p = Plots.contour(xs,ys, z,
		title="Stieltjes Transform of a White Wishart Matrix",
		size=(700,600),
		xlabel="Re(z)",
		ylabel="Im(z)",
		#zlabel="|g(z)|"
	)
	p
end

# ╔═╡ ee10d94e-9290-4d10-a794-b3e57d5c7e12
begin
	#plotlyjs()
	local z = (x,y) -> clamp(imag(g_(x,y)/π), -4,4)
	local p = Plots.surface(xs,ys, z,
		title="Stieltjes Transform of a White Wishart Matrix",
		size=(700,600),
		xlabel="Re(z)",
		ylabel="Im(z)",
		zlabel="Im g(z)"
	)
	p
end

# ╔═╡ 0a5c0410-0451-4697-bc97-51594d2a63ba
relu(x) = max(x,0.0)

# ╔═╡ 19d530ab-6943-4e42-b522-e7dd72827d6c
ρ(x,q) = (sqrt ∘ relu)((λ₊(q)- x) * (x-λ₋(q))) / (2π*q*x)

# ╔═╡ 17449982-0813-49b5-8c3f-254b20562adc
begin
	local ts = 0.0:0.01:6.0
	plot(ts, ρ.(ts,q), 
		label="",
		title = "Marcenko-Pastur Density",
		xlabel="λ",
		ylabel="ρ(λ)",
	)
end

# ╔═╡ Cell order:
# ╟─86d7fbae-2b06-4506-b635-0f3ce5cc4c4e
# ╠═b51af2cc-653e-11ee-19cc-bb1f00f82e59
# ╠═6d86c249-016e-4ac0-ae42-5ee5afa305d9
# ╠═47183ef7-a6a1-4626-8449-12950733c6cc
# ╠═84849894-e7a0-4698-a8e7-9ab08db88838
# ╠═de3d0390-844e-4053-996c-c83349b82503
# ╠═c2b60494-6de3-4e31-8dbb-911457d818e6
# ╠═e8aad398-2dff-422f-b754-2e35dea55019
# ╠═67013fbf-8b8b-4c8c-a395-9152aa760054
# ╠═62e2f4ab-9d93-49b3-a0b8-8c28156beb86
# ╠═90ff0901-8b0e-4306-b9c5-5d5bb7b2872a
# ╠═2ce15f2f-6e14-4a80-be68-acc92c7ea43f
# ╠═ee10d94e-9290-4d10-a794-b3e57d5c7e12
# ╠═0a5c0410-0451-4697-bc97-51594d2a63ba
# ╠═19d530ab-6943-4e42-b522-e7dd72827d6c
# ╠═17449982-0813-49b5-8c3f-254b20562adc
