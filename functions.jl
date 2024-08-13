using SparseArrays, StaticArrays, LinearAlgebra

function gra(cells,boxsize)
gradient = zeros(cells, cells)

for i in 1:(cells-1)
    gradient[i, i+1] = 0.5
    gradient[i+1, i] = -0.5
end
gradient[1, end] = -0.5
gradient[end, 1] = 0.5
gradient *= 1/(boxsize/cells)

gradient = SMatrix{cells, cells}(gradient)

return gradient
end

function lap(cells, boxsize)
    laplace = Matrix(I(cells))*-2
    
    for i in 1:(cells - 1)

        laplace[i, i+1] = 1
        laplace[i+1, i] = 1

    end

    laplace[end,1] = 1
    laplace[1, end] = 1

    laplace .*= (1/(boxsize/cells))^2
   # laplace = SMatrix{cells, cells}(laplace)
    return laplace

end

function gra2(cells,boxsize)
    diag1 = spdiagm(1 => fill(0.5, cells-1))
    diag2 = spdiagm(-1=> fill(-0.5,cells-1))

    corner = sparse([1, cells], [cells, 1], [-0.5, 0.5])

    mat = diag1 + diag2 + corner

    mat /= boxsize/cells

    return mat
end

function lap2(cells, boxsize)

    diag = spdiagm(0 => fill(-2, cells))
    diag1 = spdiagm(1 => fill(1, cells-1))
    diag2 = spdiagm(-1=> fill(1, cells-1))

    corner = sparse([1, cells], [cells, 1], [1, 1])

    mat = diag1 + diag2 + corner

    mat /= (boxsize/cells)^2

    return mat

end

function lap3(cells, boxsize)

#=     diag = spdiagm(0 => fill(-2, cells+2))
    diag1 = spdiagm(1 => fill(1, cells-1))
    diag2 = spdiagm(-1=> fill(1, cells-1)) =#

    mat = Tridiagonal(fill(1, cells+1), fill(-2, cells+2), fill(1, cells+1))

    mat /= (boxsize/cells)^2

    return mat

end

function circle(boxsize, particles, height)

    function circleV(x, b, y)

        top = (x-b/2)^2
        bot = (b/2)^2
        comb = top/bot
    
        right = 1-comb
    
        return sqrt(y * right)
    
    end
    
    circlePP = Vector(range(0, boxsize-boxsize/particles, particles))
    circlePN = Vector(range(boxsize/particles, boxsize, particles))
    
    circleVP = Vector([circleV(i, boxsize, height^2) for i in circlePP])
    circleVN = Vector([-circleV(i, boxsize, height^2) for i in circlePN])
    
    return append!(circlePP, circlePN), append!(circleVP, circleVN)
end    

function make_gif(P, V, frames)
   iters, particles = size(P)
   v = V[1,:]
   negativeinds = v.<0
    positiveinds = v.>=0
    anim = @animate for i in ProgressBar(1:Int(iters/frames):iters)

        Plots.scatter(P[i, negativeinds], V[i, negativeinds],
                    xlimits=(0,boxsize), 
                    ylimits=(minimum(V), maximum(V)), 
                    alpha = 0.1 + 1000/particles, 
                    markersize = 1,
                    markerstrokewidth=0,
                size=(250,250),
                # dpi = 50,
                )

        Plots.scatter!(P[i, positiveinds], V[i, positiveinds],
                    xlimits=(0,boxsize), 
                    ylimits=(minimum(V), maximum(V)), 
                    markerstrokewidth=0,
                    alpha = 0.1 + 1000/particles, 
                    markersize = 1,
                    legend=:bottomleft,
                    xlabel="Location",
                    ylabel="Speed"
                    )

    end
    gif(anim, fps = 12)
end