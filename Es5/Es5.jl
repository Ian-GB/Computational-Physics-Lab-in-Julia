### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 8b8a5edc-c863-4836-8bb3-62c97251c006
begin
	#import Pkg; Pkg.add("Plots"); Pkg.add("LaTeXStrings")
	using Random
	using Formatting
	using Plots
	theme(:ggplot2)
	using LaTeXStrings
end;

# ╔═╡ b160fd2b-4a91-4cf2-a25e-35632eb83d94
md"""
# Laboratorio di fisica computazionale
## Homework 5
### Author: Ian Gremese
#### Unresolved issues:
- Might change how IS∫ is implemented so that anything that requires a high number of repetitions will be put inside the function and be optimized as a package by the compiler. One would also want to call a "heavy" function as few times as possibile"""

# ╔═╡ 0d1f5ab3-c5f1-47e7-8ee1-c29514206574
md"""This function implements the sample mean integration algorithm."""

# ╔═╡ 95ab92f0-b18c-11ec-27bc-e1682c54b40f
function SM∫(a::Float64,b::Float64,func::Function,pts::Int64)
	
	xs = rand(Float64,pts)
	xs = xs .* (b-a) .+ a
	
	∫func = 0.::Float64
	Σᵢfᵢ² = 0.::Float64
	
	for i in xs
		fᵢ = func(i)
		∫func += fᵢ
		Σᵢfᵢ² += fᵢ^2
	end
	∫func = ∫func / pts
	Σᵢfᵢ² = Σᵢfᵢ² / pts #-> actually average of Σᵢfᵢ²
	return ∫func, Σᵢfᵢ²
	
end

# ╔═╡ 51e07c67-09d6-40eb-8342-53256b20cce6
md"""This cell contains the new compact function for ``f(x)=\exp(-x^2)``."""

# ╔═╡ 45ed1599-773e-484b-b97e-7773e62e815a
function exp⁻²(x::Float64)
	fx = Float64
	fx = exp(-x^2)
	return fx
end

# ╔═╡ 44808114-54b3-450d-a0b4-fa181dac7f88
begin
	N1 = 1000::Int64
	int = SM∫(0.,1.,exp⁻²,N1)[1]
	#dint = printfmt("{:.2f}",int)
	
	md"""In this cell we're computing the numerical estimate of ``I = \int_0^1 e^{-x^2}\mathrm d x`` with $N1 points: $(int)"""
end

# ╔═╡ 3d8cca39-350d-4e01-95f8-18f710cb1ecb
md"""In the next cell the importance sampling algorithm will be implemented."""

# ╔═╡ 2f4d934e-a05b-4cc7-a0cf-ce128347fc92
function IS∫(a::Float64,b::Float64,f::Function,pdata,pts)

	#=
	pdata should contain:
	pdata[1] = distribution function p(y)
	pdata[2] = list of random numbers according to p(y) in the interval [a,b]
	pdata[3] = ∫ₐᵇp(y)dy
	=#

	p::Function   =  pdata[1]
	prnd::Function = pdata[2]
	int::Float64  =  pdata[3]
	
	∫f = 0.::Float64
	Σᵢfᵢ² = 0.::Float64 #Where fᵢ = f(xᵢ) / p(xᵢ) * ∫ₐᵇp(x)dx
	for i = 1:1:pts
		x = prnd(a,b)
		fᵢ = f(x) / p(x)
		∫f += fᵢ
		Σᵢfᵢ² += fᵢ^2 #-> ∑ᵢ fᵢ / ∫ₐᵇp(x)dx, actually
	end
	∫f = ∫f / pts * int
	Σᵢfᵢ² = Σᵢfᵢ² / pts * int #-> <fᵢ²>, actually
	return ∫f, Σᵢfᵢ²
	
end

# ╔═╡ b0fe5f6f-57c6-4d36-9c5c-7fd4ef4f9776
md"""The next cell produces a ```pts```-component array of random numbers distributed accordingly to the exponential function ``f(x)=e^{-x}``."""

# ╔═╡ d04c8729-cd7e-4159-b4f7-40af7c4f8f8e
function exp⁻data(a::Float64,b::Float64)
	
	ν = exp(-(b-a))

	#our distribution p(x)
	function exp⁻(x::Float64)
		return exp(-x+a) / (1-ν)
	end

	#function f() that converts uniformly-distributed x into y=f(x) according to p(x)
	function e⁻rnd(a,b)
		return a - ( log(1-ν) + log1p( ν/(1-ν)  - rand(Float64) ) )::Float64
	end

	return (exp⁻, e⁻rnd, 1) #the functions have been written so that the integral is always 1
end

# ╔═╡ e8b37330-32c1-4115-9cd2-5701154d0723
IS∫(0.,1.,exp⁻²,exp⁻data(0.,1.),100000)[1]

# ╔═╡ 889a544d-a4c4-415f-bab9-cc99d56a0bc9
begin
	th = 0.746824
	pwrrng = 4:.4:8
end

# ╔═╡ 850d4796-36b5-4bbc-8e3a-51e584aa5986
md"""We are about to compute the variance (from the theoretical value) of the estimates of the integral carried out with both techniques, for a variable number ``N`` of points generated in the range $pwrrng. The variance computed as
``\sigma_{N_j} = \sum_{i=1}^{N_j} \frac{(F_n-I)^2}{N_j}``
"""

# ╔═╡ 0f74efe7-58cd-47ae-8d70-aae08b676f0f
begin
	reps = length(pwrrng)
	
	#Arrays for (Fₙ = Σᵢ fᵢ / N, Σ (Fₙ-I)² / N, Σᵢ fᵢ² / N)
	ISout = zeros(Float64,reps,3)
	SMout = zeros(Float64,reps,3)
	N = zeros(Int64,reps)

	#i = 1
	N[1] = round(10^pwrrng[1])
	#Compute estimates and add estimate to SMout[:,1] and average of square estimates to SMout[:,3]
	SMout[1,1], SMout[1,3] = SM∫(0., 1., exp⁻²,N[1])
	ISout[1,1], ISout[1,3] = IS∫(0., 1., exp⁻²,exp⁻data(0.,1.),N[1])
	#Add variance term to XXout[2]
	SMout[1,2] += ( SMout[1,1] - th )^2
	ISout[1,2] += ( ISout[1,1] - th )^2
	
	for i = 2:1:reps
		N[i] = round(10^pwrrng[i])
		#Instead of generating another N[i] values, we generate N[i]-N[i-1] and build up on the data and the "variance" computed with N[i-1]
		#In order to get some statistics, one may well make multiple attemps with the same N[i], but we are trying to highlight a trend, with no expectation to build a paper upon it.	
		#Compute estimates and add estimate to SMout[:,1] and average of square estimates to SMout[:,3]
		SMout[i,1], SMout[i,3] = SM∫(0., 1., exp⁻²,N[i]-N[i-1])
		ISout[i,1], ISout[i,3] = IS∫(0., 1., exp⁻²,exp⁻data(0.,1.),N[i]-N[i-1])
		#Process XXout[i,1] = ∑(i = N[i-1]+1, N[i]) / (N[i] - N[i-1]) to get the actual estimate
		SMout[i,1] = SMout[i,1] * (1 - N[i-1]/N[i]) + SMout[i-1,1] * N[i-1]/N[i]
		ISout[i,1] = ISout[i,1] * (1 - N[i-1]/N[i]) + ISout[i-1,1] * N[i-1]/N[i]
		#Process XXout[i,3] = ∑(i = N[i-1]+1, N[i]) / (N[i] - N[i-1]) to get the actual variance
		SMout[i,3] = SMout[i,3] * (1 - N[i-1]/N[i]) + SMout[i-1,3] * N[i-1]/N[i]
		ISout[i,3] = ISout[i,3] * (1 - N[i-1]/N[i]) + ISout[i-1,3] * N[i-1]/N[i]
		#Add variance term to XXout[2]
		SMout[i,2] += ( SMout[1,1] - th )^2
		ISout[i,2] += ( ISout[1,1] - th )^2
	end

	for i =1:1:reps
		#Rescale the variances in XXout[2]
		SMout[i,2] = sqrt(SMout[i,2] / N[i])
		ISout[i,2] = sqrt(ISout[i,2] / N[i])
		#Compute the σ = Σᵢ √(<fᵢ²> - <fᵢ>²) from <fᵢ²> (XXout[3]) and <fᵢ> (XXout[4])
		SMout[i,3] = sqrt(SMout[i,3] - (SMout[i,1])^2)
		ISout[i,3] = sqrt(ISout[i,3] - (ISout[i,1])^2)
	end
	
end

# ╔═╡ 1989fbfc-336b-4065-8e87-626e24edf958
begin
	plot(title="Uncert. on ∫₀¹ exp(-x²) dx betw. Sample Mean\n and Importance Sampling with p(x) = exp(-x).",
		yaxis=("Estimate",:log),
		xaxis=("Number of random points",:log))
	theor = fill(0.746824, length(pwrrng))
	plot!(N,ISout[:,1],label="I.S.",seriestype=:dots)
	plot!(N,SMout[:,1],label="S.M.",seriestype=:dots)
	plot!(N,theor,label="th. val.",seriestype=:path)
end

# ╔═╡ a7b7cffd-8175-4acc-bbcd-b54370e258f1
begin
	plot(title="Uncert. on ∫₀¹ exp(-x²) dx betw. Sample Mean\n and Importance Sampling with p(x) = exp(-x).",
		yaxis=("Variance on 1 attempts",:log),
		xaxis=("Number of random points",:log))
	plot!(N,ISout[:,2],label="I.S.",seriestype=:dots)
	plot!(N,SMout[:,2],label="S.M.",seriestype=:dots)
end

# ╔═╡ e1f09551-e5b9-46ab-b2a7-1c7d5bb6cb0b
begin
	plot(title="Uncert. on ∫₀¹ exp(-x²) dx betw. Sample Mean\n and Importance Sampling with p(x) = exp(-x).",
		yaxis=("Variance on 1 attempts",:log),
		xaxis=("Number of random points",:log))
	plot!(N,ISout[:,3] ./ sqrt.(N[:]),label="I.S.",seriestype=:dots)
	plot!(N,SMout[:,3] ./ sqrt.(N[:]),label="S.M.",seriestype=:dots)
end

# ╔═╡ Cell order:
# ╟─b160fd2b-4a91-4cf2-a25e-35632eb83d94
# ╠═8b8a5edc-c863-4836-8bb3-62c97251c006
# ╟─0d1f5ab3-c5f1-47e7-8ee1-c29514206574
# ╠═95ab92f0-b18c-11ec-27bc-e1682c54b40f
# ╠═51e07c67-09d6-40eb-8342-53256b20cce6
# ╠═45ed1599-773e-484b-b97e-7773e62e815a
# ╟─44808114-54b3-450d-a0b4-fa181dac7f88
# ╟─3d8cca39-350d-4e01-95f8-18f710cb1ecb
# ╟─2f4d934e-a05b-4cc7-a0cf-ce128347fc92
# ╟─b0fe5f6f-57c6-4d36-9c5c-7fd4ef4f9776
# ╟─d04c8729-cd7e-4159-b4f7-40af7c4f8f8e
# ╟─e8b37330-32c1-4115-9cd2-5701154d0723
# ╠═889a544d-a4c4-415f-bab9-cc99d56a0bc9
# ╟─850d4796-36b5-4bbc-8e3a-51e584aa5986
# ╠═0f74efe7-58cd-47ae-8d70-aae08b676f0f
# ╠═1989fbfc-336b-4065-8e87-626e24edf958
# ╠═a7b7cffd-8175-4acc-bbcd-b54370e258f1
# ╠═e1f09551-e5b9-46ab-b2a7-1c7d5bb6cb0b
