function get_acceleration(p, cells, boxsize, n0, gradient_matrix, laplace_matrix)

    n = zeros(cells)
	invdx = cells/boxsize
    dx = boxsize/cells

     for i in p
        left = Int(floor(i*invdx))
        right = left+1
        weight_left= (right * dx - i)/dx
        weight_right=(i-left*dx)/dx
        if left>=cells
            left -=cells
        end
        if right>=cells
            right -= cells
        end
        n[left+1] += weight_left
        n[right+1] += weight_right

    end

	n .*= n0 * boxsize / length(p) * invdx 

    nplus = zeros(cells+2)
    nplus[1] = n[end]
    nplus[end] = n[1]
    nplus[2:(end-1)] = n
    
    phi = (laplace_matrix \ (nplus.-n0))[2:(end-1)]

    e = vec(-gradient_matrix * phi)
    E = zeros(length(p))


     for (count,i) in enumerate(p)
        left = Int(floor(i*invdx))
        right = left+1
        weight_left= (right * dx - i)/dx
        weight_right=(i-left*dx)/dx

        if left>=cells
            left -=cells
        end

        if right>=cells
            right -= cells
        end

        a = weight_left * e[left+1]
        b = weight_right * e[right+1]

        E[count] = a + b
    end

    return -E, n, phi, e

end

function sim(cells, boxsize, particles, n0, dt, p, v, iters; all_frames = false)

    if all_frames
        P = zeros(iters, particles)
        V = zeros(iters, particles)
    end

    laplace_matrix = lap3(cells, boxsize,)
    gradient_matrix = gra2(cells, boxsize,)

    acc, n, phi, e = get_acceleration(p, cells, boxsize, n0, gradient_matrix, laplace_matrix)

    for K in ProgressBar(1:iters)

        v+= dt * acc /2
        p+= dt * v

        p[p.>boxsize] = [i%boxsize for i in p[p.>boxsize]]

        while any(p.<0)
            p[p.<0] .+= boxsize
        end

        if all_frames
            P[K, :] = p
            V[K,:] = v
        end

        acc, n, phi, e = @inline get_acceleration(p, cells, boxsize, n0, gradient_matrix, laplace_matrix)
        v = v + dt*acc / 2

    end

    #return p, v, n, phi, e, P, V
    if all_frames
        return p, v, P, V
    else
        return p, v
    end

end