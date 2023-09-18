### A Pluto.jl notebook ###
# v0.19.27

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

# ╔═╡ 937c36f6-5541-11ee-3863-d33f0332ab86
using Pkg; Pkg.activate()

# ╔═╡ f86162f7-d852-449f-b637-0f53f2cbef37
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

# ╔═╡ 7d15e3fa-cd70-449c-9110-bdf6a9595b7c
md"""
# Random Matrix Theory: Wigner Matrices and the semi-circle law

This notebook contains the solutions to some of the exercises in the second chapter of the book, but mostly some nice interactive visualizations of the Stieltjes transform.
"""

# ╔═╡ 80ee1731-4da1-47ad-8a4e-b3f239c03177
md"Comment out the following cell if you want Pluto to install the packages for you."

# ╔═╡ 90763c30-5c1e-4118-bf94-e206423169da
gr()

# ╔═╡ c528e9d5-2acb-4149-8033-6e16e27c8e2d
md"""
## Exercise 2.3.1 - Stieltjes Transform for shifted and scaled matrices
Let $A$ be a random matrix drawn from a well-behaved ensemble with Stieltjes transform $g(z)$. What are the Siteltjes transforms of the random matrices $\alpha A$ and $A + \beta I$, where $\alpha$ and $\beta$ are non-zero real numbers.
"""

# ╔═╡ 25da0e70-d4c1-4d6b-9081-011fb4de9f01
md"""
By definition of the Stieltjes transform,
```math 
g^A(z) = \frac{1}{N} \text{tr}((zI - A)^{-1} ),
```
and $(zI - \alpha A) = \alpha (\frac{z}{\alpha} I - A)$
hence
```math
g^{\alpha A}(z) = \frac{1}{N} \frac{1}{\alpha} \text{tr} ((\frac{z}{\alpha}I - A)^{-1}) = \frac{1}{\alpha} g^A(\frac{z}{\alpha})
```
Similarly, 
```math
g^{A + \beta I}(z) = \frac{1}{N} \text{tr} (((z-\beta)I - A)^{-1}) = g(z-\beta)
```
"""

# ╔═╡ a6e30f43-809d-4996-925c-22a5d9251201
stieltjes_transform(λs,z) = mean(1/(z-λ) for λ in λs)

# ╔═╡ 3ed04034-42b1-4b9c-bbb5-8d06a9dc3996
stieltjes_transform(A::AbstractMatrix,z) = stieltjes_transform(eigvals(A), z)

# ╔═╡ 1738150e-a074-4c0a-82b0-2759d39858d5
function random_wigner(n, σ=1.0)
	H = randn(n,n)
	return Symmetric((H+H') * σ / sqrt(2*n))
end

# ╔═╡ ae5ba9db-5277-4600-8cf4-8cf40ddb2688
begin
	xmin = -4.0
	xmax = 4.0
	ymin = -4.0
	ymax = 4.0
	step = 0.01
	xs = xmin:step:xmax
	ys = ymin:step:ymax
end

# ╔═╡ cb27a1e8-c440-4f48-adba-cd873fa23d3b
md"""
You can move the slider below to change the size of the matrix used for the plot.
"""

# ╔═╡ a28aa970-8441-4b68-9d95-d1c804f7d4c4
@bind n Slider(1:100; default=5, show_value=true)

# ╔═╡ 3d552c8d-88af-4e60-b82a-f55492d0957c
begin
	plotlyjs() # switch to plotlyjs backend for this plot
	local A = random_wigner(n)
	local λs = eigvals(A)
	local z = (x,y) -> min(4, abs(stieltjes_transform(λs, x+y*im)))
	p = Plots.surface(xs,ys, z,
		title="Stieltjes Transform of a Wigner Matrix (N=$n)",
		size=(700,600),
		xlabel="Re(z)",
		ylabel="Im(z)",
		zlabel="|g(z)|"
	)
	local η_1 = fill(-1/n, length(xs))
	local η_2 = fill(-1/sqrt(n), length(xs))
	#scatter!(p, xs, η_1, z.(xs,η_1), label=L"x-i/N")
	#scatter!(p, xs, η_2, z.(xs,η_2), label=L"x-i/\sqrt{N}")
	p
end

# ╔═╡ c23aec92-41a6-4071-924b-1215d05185f2
g_(z::Complex,σ=1) = (z - z*sqrt(1 -4(σ^2)/z^2))/(2σ^2)

# ╔═╡ 488e41f9-8ceb-41da-a709-075cbb20e162
g_(x,y) = g_(x+y*im)

# ╔═╡ 770b9719-ad78-460a-8098-9b9ef494c6e4
begin
	plotlyjs()
	local z = (x,y) -> min(4, abs(g_(x,y)))
	local p = Plots.surface(xs,ys, z,
		title="Stieltjes Transform of a Wigner Matrix (Analytical limit)",
		size=(700,600),
		xlabel="Re(z)",
		ylabel="Im(z)",
		zlabel="|g(z)|"
	)
	p
end

# ╔═╡ aac7f646-d0d7-43d1-a4d2-12669d39ddd8
md"""
Note that it is different than in the video! Not sure if it was a typo or just complex functions being finicky.
"""

# ╔═╡ 9dc436c2-f7f0-42a3-96b7-adb4b4bad369
md"""
Here's the plot of the imaginary part. We can clearly see how it approaches the semi-circle density as we approach from the real axis from below.
"""

# ╔═╡ 9c0ebc3a-225c-454d-ac18-33810e6f33b5
begin
	plotlyjs()
	local z = (x,y) -> min(4, imag(g_(x,y)))
	local p = Plots.surface(xs,ys, z,
		title="Stieltjes Transform of a Wigner Matrix (Analytical limit)",
		size=(700,600),
		xlabel="Re(z)",
		ylabel="Im(z)",
		zlabel="Im g(z)"
	)
	p
end

# ╔═╡ d0dbbe01-af4e-44b8-8c64-ecde5b5e4242
md"""
## Exercise 2.3.2 Finite ``N`` approximation and small imaginary part
``\text{Im} g_N(x-i\eta)/\pi`` is a good approximation of ``\rho(x)`` for small positive ``\eta`` where ``g_N(z)`` is the sample Stieltjes transform.

Numerically generate a Wigner matrix of size ``N`` and ``\sigma^2=1``.

a) For ``\eta \in \{ 1/N, 1/\sqrt{N}, 1 \}``, plot $\text{Im} g_N(x-i\eta)/\pi$ and the theoretical $\rho(x)$ on the same plot.
"""

# ╔═╡ 2a24ce45-563e-4a42-81ee-e3d873860ffd
ρ(x, σ=1) = abs(x) > 2σ ? 0.0 : sqrt(4σ^2 - x^2)/(2π*σ^2)

# ╔═╡ 31051193-340d-4893-9850-04aebdf2d0cc
ρ_approx(λ) = (x,y) -> imag(stieltjes_transform(λ, x-y*im)) / pi

# ╔═╡ f1c34391-de42-4d7d-bbcf-006b00a377e3
begin
	gr()
	local n = 400
	local A = random_wigner(n)
	local λs = eigvals(A)
	local g = ρ_approx(λs)
	local p = plot(xs, ρ, label="Analytical density")
	plot!(p, xs, g.(xs, 1.0), label=L"\eta=1")
	plot!(p, xs, g.(xs, 1/sqrt(n)), label=L"\eta=1/\sqrt{N}")
	#plot!(p, xs, g.(xs, 1/n), label=L"\eta=1/N")
	p
end

# ╔═╡ 5b19482a-883c-4183-99f3-18229b10fc1e
md"""
b) Compute the error as a function of ``\eta`` where the error is 
```math
(\rho(x)-\text{Im}g_N(x-i\eta)/\pi)^2
```
integrated over ``[-3,3]`` (use equidistant point with $0.01$ spacing to approximate the integral). Plot this error for $\eta$ between $1/N$ and $1$. You should see that $1/\sqrt{N}$ is close to the minimum of this function.
"""

# ╔═╡ afb6fd9d-3e91-4c58-a96f-a6dbaa58e532
xs_err = -3:0.01:3

# ╔═╡ f39bf84c-b21a-475a-ae34-4379fe07f4c4
err_fun(g,x,η) = (ρ(x)-g(x,η))^2

# ╔═╡ 890230ae-216f-4d2d-b895-b7fbb770dcba
err(g,η) = sum(err_fun(g,x,η) for x in xs_err)

# ╔═╡ c3588633-07fc-4cf1-a08d-f806c7817701
begin
	local n = 400
	local A = random_wigner(n)
	local λs = eigvals(A)
	local g = ρ_approx(λs)
	ηs = (1/n):0.01:1
	plot(ηs, err.(g, ηs), label="")
	vline!([1/sqrt(n)], label=L"1/\sqrt{N}")
end

# ╔═╡ 342379ba-be1f-4362-9cf3-b6b8edeaaa31
md"""
## Ginibre Ensemble
Let $X \in \mathbb{R}^{N\times N}$ be a random matrix with $X_{ij} \sim N(0,\sigma^2/N)$, then the eigenvalue density of $X$ follows the Girko circular law. Eigenvalues are uniformly distributed in a circle centered at $0$ and of radius $\sigma^2$
"""

# ╔═╡ cf025540-4644-4a2d-8868-bfcc63614d27
begin
	local n = 1000
	local A = randn(n,n) / sqrt(n)
	local λs = eigvals(A)
	scatter(λs, label="", title="Girko circular law")
	local t = 0.0:0.01:2π
	plot!(cos.(t), sin.(t), label="", linestyle=:dash)
end

# ╔═╡ ea1a543c-7896-46c0-8f47-1cb35bb3f69f
begin
	local n = 1000
	local A = randn(n,n) / sqrt(2n)
	local B = randn(n,n) / sqrt(2n)
	local λs = eigvals(A+im*B)
	scatter(λs, label="", title="Girko circular law (Complex Gaussian)")
	local t = 0.0:0.01:2π
	plot!(cos.(t), sin.(t), label="", linestyle=:dash)
end

# ╔═╡ 1338402d-0905-46a2-95a8-42c3f50d5918
# ╠═╡ disabled = true
#=╠═╡
md"""
## Exercise 2.3.3 From the moments to the density
"""
  ╠═╡ =#

# ╔═╡ 70ed779e-417f-48b4-899b-4df3f59db796
# ╠═╡ disabled = true
#=╠═╡
md"""
## Exercise 3.1.1 Quaternionic matrices of size one
"""
  ╠═╡ =#

# ╔═╡ ce009503-175b-4c74-94e6-e8036f9ff6b3
# ╠═╡ disabled = true
#=╠═╡
md"""
## Exercise 3.1.2 Three quarter-circle laws
"""
  ╠═╡ =#

# ╔═╡ 401a26a1-8363-45c6-9a87-58ee1bf4ff55
# ╠═╡ disabled = true
#=╠═╡
md"""
## Exercise 3.2.1 Non-crossing pair partitions of eight elements
"""
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─7d15e3fa-cd70-449c-9110-bdf6a9595b7c
# ╟─80ee1731-4da1-47ad-8a4e-b3f239c03177
# ╠═937c36f6-5541-11ee-3863-d33f0332ab86
# ╠═f86162f7-d852-449f-b637-0f53f2cbef37
# ╠═90763c30-5c1e-4118-bf94-e206423169da
# ╟─c528e9d5-2acb-4149-8033-6e16e27c8e2d
# ╟─25da0e70-d4c1-4d6b-9081-011fb4de9f01
# ╠═a6e30f43-809d-4996-925c-22a5d9251201
# ╠═3ed04034-42b1-4b9c-bbb5-8d06a9dc3996
# ╠═1738150e-a074-4c0a-82b0-2759d39858d5
# ╠═ae5ba9db-5277-4600-8cf4-8cf40ddb2688
# ╟─cb27a1e8-c440-4f48-adba-cd873fa23d3b
# ╠═a28aa970-8441-4b68-9d95-d1c804f7d4c4
# ╠═3d552c8d-88af-4e60-b82a-f55492d0957c
# ╠═c23aec92-41a6-4071-924b-1215d05185f2
# ╠═488e41f9-8ceb-41da-a709-075cbb20e162
# ╠═770b9719-ad78-460a-8098-9b9ef494c6e4
# ╟─aac7f646-d0d7-43d1-a4d2-12669d39ddd8
# ╟─9dc436c2-f7f0-42a3-96b7-adb4b4bad369
# ╠═9c0ebc3a-225c-454d-ac18-33810e6f33b5
# ╟─d0dbbe01-af4e-44b8-8c64-ecde5b5e4242
# ╠═2a24ce45-563e-4a42-81ee-e3d873860ffd
# ╠═31051193-340d-4893-9850-04aebdf2d0cc
# ╠═f1c34391-de42-4d7d-bbcf-006b00a377e3
# ╟─5b19482a-883c-4183-99f3-18229b10fc1e
# ╠═afb6fd9d-3e91-4c58-a96f-a6dbaa58e532
# ╠═f39bf84c-b21a-475a-ae34-4379fe07f4c4
# ╠═890230ae-216f-4d2d-b895-b7fbb770dcba
# ╠═c3588633-07fc-4cf1-a08d-f806c7817701
# ╟─342379ba-be1f-4362-9cf3-b6b8edeaaa31
# ╠═cf025540-4644-4a2d-8868-bfcc63614d27
# ╠═ea1a543c-7896-46c0-8f47-1cb35bb3f69f
# ╟─1338402d-0905-46a2-95a8-42c3f50d5918
# ╟─70ed779e-417f-48b4-899b-4df3f59db796
# ╟─ce009503-175b-4c74-94e6-e8036f9ff6b3
# ╟─401a26a1-8363-45c6-9a87-58ee1bf4ff55
