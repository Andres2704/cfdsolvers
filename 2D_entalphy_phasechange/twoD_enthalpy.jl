module twoD_enthalpy

using LinearAlgebra

export solve_enthalpy

function initialize(domain, mesh, ics, mat_prop)
    # Domain and mesh definition
    L = domain.L 
    hi = domain.hi 
    ha = domain.ha 
    nx = mesh.NX 
    nz = mesh.NZ

    # Initialization of the matrices
    H = zeros(Float64, nz, nx)
    T = zeros(Float64, nz, nx)
    rho = zeros(Float64, nz, nx)
    S = zeros(Float64, nz, nx)
    lambda = zeros(Float64, nz, nx)  
    phi = zeros(Float64, nz, nx)   

    # Entalphy of domain e
    if ics.Te_0 < mat_prop.T_m
        He = mat_prop.cp_i * (ics.Te_0 - mat_prop.T_m)
    else
        He = mat_prop.cp_w * (ics.Te_0 - mat_prop.T_m) + mat_prop.Lm
    end

    x = collect(LinRange(-0.5*L, 0.5*L, nx))                # x domain
    z = collect(LinRange(-ha, hi, nz))                      # z domain 
    X = x' .* ones(length(z))                               # X - grid 
    Z = ones(length(x))'.* z                                # Z - grid 
    
    # Initialization of domain a and e 
    for i in eachindex(z)
        if z[i] < 0.0
            S[i,:] .= ics.S
            rho[i,:] .= mat_prop.rhoa 
            T[i, :] .= ics.Ta_0
            H[i, :] .= mat_prop.cp_a*(ics.Ta_0 - mat_prop.T_m) + mat_prop.Lm
            lambda[i,:] .= mat_prop.l_a
            phi[i, :] .= calculate_phi(z[i], H[i, 1], mat_prop.Lm) 
        else
            rho[i,:] .= mat_prop.rhoe 
            T[i, :] .= ics.Te_0
            H[i, :] .= He
            lambda[i,:] .= lambda_e(He, mat_prop)
            phi[i, :] .= calculate_phi(z[i], He, mat_prop.Lm)
        end
    end

    return x, z, X, Z, H, T, rho, S, lambda, phi
end

function lambda_e(He, mat_prop)
    if He<0
        return mat_prop.l_i
    elseif He > mat_prop.Lm
        return mat_prop.l_w
    else 
        phi = He/mat_prop.Lm
        return (phi/mat_prop.l_w + (1-phi)/mat_prop.l_i)^(-1)
    end    
end

function calculate_temperature(z, H, mat_prop)
    if z < 0.0
        T = (H -mat_prop.Lm)/mat_prop.cp_a + mat_prop.T_m
    else
        if (H<0.0)
            T = H/mat_prop.cp_i + mat_prop.T_m 
        elseif (H <= mat_prop.Lm) && (H>=0.0)
            T = mat_prop.T_m
        else
            T = (H - mat_prop.Lm)/mat_prop.cp_w + mat_prop.T_m 
        end
    end
    return T 
end

function boundary_conditions(z, H, T, rho, S, lambda, dx, dz, dt, nx, nz, bcs)
    # ================================== Boundary conditions ================================
    for j in 2:nx-1
        # Top Boundary
        H[nz,j] = H[nz,j] -T[nz,j]*(dt/rho[nz,j]) * ( (lambda[nz,j+1] + 2*lambda[nz,j] + lambda[nz,j-1])/(2*dx^2) +
                                               (3*lambda[nz,j] + lambda[nz-1,j])/(2*dz^2) + bcs.htc/dz) + 
                    T[nz,j+1]*(dt/rho[nz,j]) * (lambda[nz,j+1] + lambda[nz,j])/(2*dx^2) + 
                    T[nz,j-1]*(dt/rho[nz,j]) * (lambda[nz,j] + lambda[nz,j-1])/(2*dx^2) +
                    T[nz-1,j]*(dt/rho[nz,j]) * (3*lambda[nz,j] + lambda[nz-1,j])/(2*dz^2) + bcs.htc*bcs.Trec/dz + S[nz,j]
        
        # Bottom bord 
        H[1,j] = H[1,j] -T[1,j]*(dt/rho[1,j]) * ( (lambda[1,j+1] + 2*lambda[1,j] + lambda[1,j-1])/(2*dx^2) +
                                                  (lambda[2,j] + 3*lambda[1,j])/(2*dz^2) ) + 
                    T[1,j+1]*(dt/rho[1,j]) * (lambda[1,j+1] + lambda[1,j])/(2*dx^2) + 
                    T[1,j-1]*(dt/rho[1,j]) * (lambda[1,j] + lambda[1,j-1])/(2*dx^2) +
                    T[2,j  ]*(dt/rho[1,j]) * (lambda[2,j] + 3*lambda[1,j])/(2*dz^2) + S[1,j]
    end 

    for i in 2:nz-1
        # Left and Right bord 
        H[i,1] = H[i,1] -T[i,1]*(dt/rho[i,1]) * ( (lambda[i,2] + 3*lambda[i,1])/(2*dx^2) +
                                                  (lambda[i+1,1] + 2*lambda[i,1] + lambda[i-1,1])/(2*dz^2) ) + 
                        T[i,2]*(dt/rho[i,1]) * (lambda[i,2] + 3*lambda[i,1])/(2*dx^2) + 
                        T[i+1,1]*(dt/rho[i,1]) * (lambda[i+1,1] + lambda[i,1])/(2*dz^2) + 
                        T[i-1,1]*(dt/rho[i,1]) * (lambda[i,1] + lambda[i-1,1])/(2*dz^2) + S[i,1]

        H[i,nx] = H[i,nx] -T[i,nx]*(dt/rho[i,nx]) * ( (3*lambda[i,nx] + lambda[i,nx-1])/(2*dx^2) +
                                                (lambda[i+1,nx] + 2*lambda[i,nx] + lambda[i-1,nx])/(2*dz^2) ) + 
                        T[i,nx-1]*(dt/rho[i,nx]) * (3*lambda[i,nx] + lambda[i,nx-1])/(2*dx^2) +
                        T[i+1,nx]*(dt/rho[i,nx]) * (lambda[i+1,nx] + lambda[i,nx])/(2*dz^2) + 
                        T[i-1,nx]*(dt/rho[i,nx]) * (lambda[i,nx] + lambda[i-1,nx])/(2*dz^2) + S[i,nx]
    end     

    return H     
end

function find_front(z, H, nx, mat_prop)
    front_eight = Float64[]
    for i in eachindex(z)
        if (z[i]>=0.0)
            phii = calculate_phi(z[i], H[i, Int64(nx/2)], mat_prop.Lm)
            if (phii >= 0.0) && (phii <= 1)
                push!(front_eight, z[i])
            end
        end 
    end
    return minimum(front_eight)
end

function calculate_phi(z, He, Lm)
    if He <= 0.0
            phi = 0.0
    elseif  (He <= Lm) && (He >= 0.0)
            phi = He/Lm
    else
            phi = 1
    end
    return phi 
end

function propagate_solution!(z, H, T, rho, S, lambda, phi, dx, dz, nx, nz, dt, mat_prop, bcs)
    H = boundary_conditions(z, H, T, rho, S, lambda, dx, dz,  dt, nx, nz, bcs)    
    for i in 2:(nz-1)
        for j in 2:(nx-1)
            H[i,j] = H[i,j] -T[i,j]*(dt/rho[i,j]) * ( (lambda[i,j+1] + 2*lambda[i,j] + lambda[i,j-1])/(2*dx^2) +
                                                      (lambda[i+1,j] + 2*lambda[i,j] + lambda[i-1,j])/(2*dz^2) ) + 
                    T[i,j+1]*(dt/rho[i,j]) * (lambda[i,j+1] + lambda[i,j])/(2*dx^2) + 
                    T[i,j-1]*(dt/rho[i,j]) * (lambda[i,j] + lambda[i,j-1])/(2*dx^2) +
                    T[i+1,j]*(dt/rho[i,j]) * (lambda[i+1,j] + lambda[i,j])/(2*dz^2) + 
                    T[i-1,j]*(dt/rho[i,j]) * (lambda[i,j] + lambda[i-1,j])/(2*dz^2) + S[i,j]
        end
    end

    for i in 1:nz
        for j in 1:nx
            # Updating the temperature based on T-h relation 
            T[i,j] = calculate_temperature(z[i], H[i,j], mat_prop)
            phi[i,j] = calculate_phi(z[i], H[i,j], mat_prop.Lm)

            # Updating the thermal conductivity based on enthalpy
            if z[i] < 0.0
                lambda[i,j] = mat_prop.l_a 
            else
                lambda[i,j] = lambda_e(H[i,j], mat_prop)
            end 
        end
    end

    return H, T, rho, S, lambda, phi 
end

function propagate_solution_newton!(z, H, T, rho, S, lambda, phi, dx, dz, nx, nz, dt, mat_prop, bcs;
                                    max_iter=20, tol=1e-6)

    H_old = copy(H)  # Valeur connue à tⁿ⁻¹

    for iter = 1:max_iter
        # === 1. Mise à jour de T, λ, ϕ à partir de H actuel ===
        for i in 1:nz
            for j in 1:nx
                T[i,j] = calculate_temperature(z[i], H[i,j], mat_prop)
                phi[i,j] = calculate_phi(z[i], H[i,j], mat_prop.Lm)

                if z[i] < 0.0
                    lambda[i,j] = mat_prop.l_a
                else
                    lambda[i,j] = lambda_e(H[i,j], mat_prop)
                end
            end
        end

        # === 2. Appliquer les conditions aux limites (bords) ===
        H = boundary_conditions(z, H, T, rho, S, lambda, dx, dz,  dt, nx, nz, bcs)    
        for i in 2:(nz-1)
            for j in 2:(nx-1)
                H[i,j] = H_old[i,j] -T[i,j]*(dt/rho[i,j]) * ( (lambda[i,j+1] + 2*lambda[i,j] + lambda[i,j-1])/(2*dx^2) +
                                                        (lambda[i+1,j] + 2*lambda[i,j] + lambda[i-1,j])/(2*dz^2) ) + 
                        T[i,j+1]*(dt/rho[i,j]) * (lambda[i,j+1] + lambda[i,j])/(2*dx^2) + 
                        T[i,j-1]*(dt/rho[i,j]) * (lambda[i,j] + lambda[i,j-1])/(2*dx^2) +
                        T[i+1,j]*(dt/rho[i,j]) * (lambda[i+1,j] + lambda[i,j])/(2*dz^2) + 
                        T[i-1,j]*(dt/rho[i,j]) * (lambda[i,j] + lambda[i-1,j])/(2*dz^2) + S[i,j]
            end
        end

        # === 4. Approximation de la jacobienne numériquement ===
        δH = zeros(nz, nx)
        ε = 1e-6
        H_pert = copy(H)
        H_pert .+= ε

        # === 3. Calcul du résidu (matrice R même taille que H) ===
        R = zeros(nz, nx)
        for i in 2:(nz-1)
            for j in 2:(nx-1)
                R[i,j] = rho[i,j]*(H[i,j] - H_old[i,j])/dt +
                         T[i,j]*((lambda[i,j+1] + 2*lambda[i,j] + lambda[i,j-1])/(2*dx^2) +
                                 (lambda[i+1,j] + 2*lambda[i,j] + lambda[i-1,j])/(2*dz^2)) -
                         T[i,j+1]*(lambda[i,j+1] + lambda[i,j])/(2*dx^2) -
                         T[i,j-1]*(lambda[i,j]   + lambda[i,j-1])/(2*dx^2) -
                         T[i+1,j]*(lambda[i+1,j] + lambda[i,j])/(2*dz^2) -
                         T[i-1,j]*(lambda[i,j]   + lambda[i-1,j])/(2*dz^2) - S[i,j]

                # Recalcule T et λ localement
                T_pert = calculate_temperature(z[i], H_pert[i,j], mat_prop)
                lambda_pert = (z[i] < 0.0) ? mat_prop.l_a : lambda_e(H_pert[i,j], mat_prop)

                # Estimation résidu perturbé
                R_pert = rho[i,j]*(H_pert[i,j] - H_old[i,j])/dt +
                         T_pert*((lambda[i,j+1] + 2*lambda_pert + lambda[i,j-1])/(2*dx^2) +
                                 (lambda[i+1,j] + 2*lambda_pert + lambda[i-1,j])/(2*dz^2)) -
                         T[i,j+1]*(lambda[i,j+1] + lambda_pert)/(2*dx^2) -
                         T[i,j-1]*(lambda_pert + lambda[i,j-1])/(2*dx^2) -
                         T[i+1,j]*(lambda[i+1,j] + lambda_pert)/(2*dz^2) -
                         T[i-1,j]*(lambda_pert + lambda[i-1,j])/(2*dz^2) - S[i,j]

                δH[i,j] = -R[i,j] / ((R_pert - R[i,j]) / ε)  # Estimation de dF/dH⁻¹ * -F
            end
        end

        res_norm = norm(R)
        # if res_norm < tol
        #     println("Newton convergé en $iter itérations avec ||F|| = $res_norm")
        #     break
        # elseif iter == max_iter
        #     println("⚠️ Newton n’a pas convergé après $max_iter itérations (||F|| = $res_norm)")
        # end

        # === 5. Mise à jour ===
        H .+= 1.879 * δH

    end 

    return H, T, rho, S, lambda, phi
end


function solve_enthalpy(domain, mesh, mat_prop, bcs, ics, solver)
    # Domain and mesh definition
    L = domain.L 
    hi = domain.hi 
    ha = domain.ha 
    nx = mesh.NX 
    nz = mesh.NZ
    dx, dz = L / nx, (hi + ha) / nz                       # mesh step-size
    
    front = Float64[0.0]

    x, z, X, Z, H, T, rho, S, lambda, phi  = initialize(domain, mesh, ics, mat_prop)

    t = 0.0
    dt = 5e-3
    while (t < solver.END_TIME)
       
       H, T, rho, S, lambda, phi = propagate_solution_newton!(z, H, T, rho, S, lambda, phi, dx, dz, nx, nz, dt, mat_prop, bcs) 
       
       push!(front, find_front(z, H, nz, mat_prop))

       t = t + dt 

    end

    return x, z, X, Z, H, T, rho, S, lambda, front, phi
end

end