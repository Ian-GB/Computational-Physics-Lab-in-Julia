### A Pluto.jl notebook ###
# v0.18.2

using Markdown
using InteractiveUtils

# ╔═╡ 8b8a5edc-c863-4836-8bb3-62c97251c006
begin
	using PlutoUI
	using Random
	using Formatting
	using Plots
end

# ╔═╡ 0d1f5ab3-c5f1-47e7-8ee1-c29514206574
md"""This function implements the sample mean integration algorithm."""

# ╔═╡ 95ab92f0-b18c-11ec-27bc-e1682c54b40f
function SM∫(a::Float64,b::Float64,func::Function,pts::Int64)
	
	xs = rand(Float64,pts)
	xs = xs .* (b-a) .+ a
	∫func = Float64

	∫func = 0
	for i in xs
		∫func += func(i)
	end
	∫func = ∫func / pts

	return ∫func
	
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
	N1 = 100
	int = SM∫(0.,1.,exp⁻²,100)
	#dint = printfmt("{:.2f}",int)
	
	md"""In this cell we're computing the numerical estimate of ``I = \int_0^1 e^{-x^2}\mathrm d x`` with $N1 points: $(int)"""
end

# ╔═╡ 3d8cca39-350d-4e01-95f8-18f710cb1ecb
md"""In the next cell the importance sampling algorithm will be implemented."""

# ╔═╡ 2f4d934e-a05b-4cc7-a0cf-ce128347fc92
function IS∫(a::Float64,b::Float64,f::Function,pdata,pts::Int64)

	#=
	pdata should contain:
	pdata[1] = distribution function p(y)
	pdata[2] = function f() that converts uniformly-distributed x into y=f(x) distributed according to p(y) in the interval [a,b]
	pdata[3] = ∫ₐᵇp(y)dy
	=#

	p::Function   =  pdata[1]
	prnd::Function = pdata[2]
	int::Float64  =  pdata[3]
	
	xs = prnd(a,b,pts) #Generates random xs according to distribution p(x)
	
	∫f = 0.::Float64
	for i in xs
		∫f += f(i) / p(i)
	end
	∫f = ∫f / pts
	∫f = ∫f * int

	return ∫f
	
end

# ╔═╡ b0fe5f6f-57c6-4d36-9c5c-7fd4ef4f9776
md"""The next cell produces a ```pts```-component array of random numbers distributed accordingly to the exponential function ``f(x)=e^{-x}``."""

# ╔═╡ d04c8729-cd7e-4159-b4f7-40af7c4f8f8e
function expdata(a,b)
	
	function exp⁻(x)
		return exp(-x)
	end

	function e⁻rnd(a::Float64,b::Float64,pts::Int64)
		x = zeros(Float64,pts)
		ν = exp(-(b-a))
		for i = 1:1:pts
			x[i] = a - ( log(1-ν) + log1p( ν/(1-ν)  - rand(Float64) ) )
		end
		return x
	end

	return (exp⁻, e⁻rnd, 1) #the functions have been written so that the integral is always 1
end

# ╔═╡ e8b37330-32c1-4115-9cd2-5701154d0723
IS∫(0.,1.,exp⁻²,expdata(0.,1.),1000)

# ╔═╡ cfe4cea2-e822-4d5f-9e22-336f9612803b
begin
	df = expdata(0.,1.)
	df[2](0.,1.,100)
	plot(x,seriestype = :histogram)
end

# ╔═╡ Cell order:
# ╠═8b8a5edc-c863-4836-8bb3-62c97251c006
# ╟─0d1f5ab3-c5f1-47e7-8ee1-c29514206574
# ╠═95ab92f0-b18c-11ec-27bc-e1682c54b40f
# ╠═51e07c67-09d6-40eb-8342-53256b20cce6
# ╠═45ed1599-773e-484b-b97e-7773e62e815a
# ╠═44808114-54b3-450d-a0b4-fa181dac7f88
# ╠═3d8cca39-350d-4e01-95f8-18f710cb1ecb
# ╠═2f4d934e-a05b-4cc7-a0cf-ce128347fc92
# ╟─b0fe5f6f-57c6-4d36-9c5c-7fd4ef4f9776
# ╠═d04c8729-cd7e-4159-b4f7-40af7c4f8f8e
# ╠═e8b37330-32c1-4115-9cd2-5701154d0723
# ╠═cfe4cea2-e822-4d5f-9e22-336f9612803b
