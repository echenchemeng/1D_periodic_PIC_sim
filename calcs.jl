using LinearAlgebra
using Plots, Statistics, Random
using NearestNeighbors
using ProgressBars, BenchmarkTools

include("functions.jl")
include("acceleration.jl")

Random.seed!(4)

cells= 1000
boxsize= 50
particles = 50000
dt = 1
iters = 50
n0 = 1

p = (zeros(particles))
v = (zeros(particles))


#= p = rand(particles)*boxsize =#

p[1:(Int(particles/2))] = range(0, boxsize - boxsize/cells, Int(particles/2))
p[(Int(particles/2)+1):end] = range(0, boxsize - boxsize/cells, Int(particles/2))


#= v[(Int(particles/2)+1):end] = range(-1, 1, Int(particles/2)) .-3
v[1:(Int(particles/2))] = range(-1, 1, Int(particles/2)) .+3 =#

v[(Int(particles/2)+1):end] .= -3
v[1:(Int(particles/2))] .= +3


#=  v[(Int(particles/2)+1):end] = randn(Int(particles/2)) .-3
 v[1:(Int(particles/2))] = randn(Int(particles/2)) .+3 =#


v .*= (1 .+ 0.2 .* map(sin, 2*pi * p / boxsize))

negativeinds = v.<0
positiveinds = v.>=0


fig = Plots.scatter(p[negativeinds],v[negativeinds],
                xlimits=(0,boxsize), alpha = 0.05, markersize = 1)
Plots.scatter!(p[positiveinds],v[positiveinds],
                xlimits=(0,boxsize), alpha = 0.05, markersize = 1)



p, v,P,V = @inline sim(cells, boxsize, particles, n0, dt, p, v, iters, all_frames=true)


Plots.scatter!(p[negativeinds],v[negativeinds],
                xlimits=(0,boxsize), alpha = 0.2, markersize = 2)
Plots.scatter!(p[positiveinds],v[positiveinds],
                xlimits=(0,boxsize), alpha = 0.1, markersize = 2)
plot!(size=(400,400), legend=false,
        xlabel="Location",
        ylabel="Speed"
        
        )


Plots.display(fig)

make_gif(P,V, 10)
