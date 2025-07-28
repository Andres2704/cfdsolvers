module cavity_functions_rhiechow

export solve_cavity

function initialize_variables(nx, ny)
    #= 
    Function used to initialize variables for a given grid.

    Inputs: 
    - nx: Number of internal elements in x 
    - ny: Number of internal elements in y 

    Outputs: 
    - u: u-velocity centered on the face (left face)
    - v: v-velocity centered on the face (bottom face)
    - p: Pressure centered in the cell
    - div: Divergence centered in the cell
    - u_old: Velocity at previous time step
    - v_old: Velocity at previous time step
    =# 
    u = zeros(Float64, ny+2, nx+2)
    v = zeros(Float64, ny+2, nx+2)
    p = zeros(Float64, ny+2, nx+2)
    div = zeros(Float64, ny+2, nx+2)
    u_old = zeros(Float64, ny+2, nx+2)
    v_old = zeros(Float64, ny+2, nx+2)

    return u, v, p, div, u_old, v_old
end

function boundary_conditions!(u_old, v_old, ut, ub, ul, ur, vt, vb, vl, vr, nx, ny)
    #=
    Function to apply boundary conditions based on given wall velocities

    Inputs:
    - u_old: Known u-velocity centered on the face 
    - v_old: Known v-velocity centered on the face
    - ut, vt: u and v velocity at top face
    - ub, vb: u and v velocity at bottom face
    - ul, vl: u and v velocity at left face 
    - ur, vr: u and v velocity at right face
    - nx: Number of internal elements in x 
    - ny: Number of internal elements in y 

    Outputs:
    - u_old: u-velocity with proper boundary conditions applied
    - v_old: v-velocity with proper boundary conditions applied
    =#

    # Applying BC's on u-velocity component
    u_old[:, 2] .= ul
    u_old[:, nx+2] .= ur
    u_old[ny+2, :] .= 2*ut .- u_old[ny+1, :]
    u_old[1, :] .= 2*ub .- u_old[2, :]

    # Applying BC's on v-velocity component 
    v_old[ny+2, :] .= vt
    v_old[2, :] .= vb
    v_old[:, 1] .= 2*vl .- v_old[:, 2]
    v_old[:, nx+2] .= 2*vr .- v_old[:, nx+1]

    return u_old, v_old, u_old, v_old
end

function stability_condition(u, v, dx::Real, dy::Real; CFL::Real = 0.5)
    #=
    Function to calculate the time step that respects the stability condition

    Inputs:
    - u: u-velocities 
    - v: v-velocities 
    - dx: Spatial step in x coordinate
    - dy: Spatial step in y coordinate 
    - CFL: Desired CFL condition 

    Outputs: 
    - dt: Time step based on advection equation stability condition
    =# 

    max_vel = maximum(sqrt.(u.^2 + v.^2))
    dt = CFL * minimum([dx, dy]) / (max_vel + 1e-10)
    return dt 
end

function intermediate_velocity!(u, v, u_old, v_old, rho, mu, nx, ny, dx, dy, dt)
    #=
    Function to solve the conservation equation for the intermediate step (see lecture notes)

    Inputs:
    - u: u-velocity centered on the face - unknown 
    - v: v-velocity centered on the face - unknown
    - u_old: Known u-velocity centered on the face 
    - v_old: Known v-velocity centered on the face
    - rho: Fluid density [SI]
    - mu: Dynamic viscosity of fluid [SI]
    - nx: Number of internal elements in x 
    - ny: Number of internal elements in y 
    - dx: Step in x direction
    - dy: Step in y direction 
    - dt: Time step 

    Outputs:
    - u: u-velocity at intermediate step - Solution
    - v: v-velocity at intermediate step - Solution
    =#

    # Solution for u-velocity 
    for i in 2:nx+1 # Iteration in x direction
        for j in 2:ny+1 # Iteration in y direction
            ue = 0.5 * (u_old[j, i]   + u_old[j, i+1])
            uw = 0.5 * (u_old[j, i-1] + u_old[j, i])
            un = 0.5 * (u_old[j+1, i] + u_old[j, i])
            us = 0.5 * (u_old[j-1, i] + u_old[j, i])
            vn = 0.5 * (v_old[j+1, i-1] + v_old[j+1, i])
            vs = 0.5 * (v_old[j, i-1]   + v_old[j, i])

            # Convective term 
            convection = -(rho * ue^2 - rho * uw^2) / dx - (rho * vn * un - rho * vs * us) / dy

            # Diffusive term 
            diffusion = mu * (u_old[j, i-1] - 2*u_old[j, i] + u_old[j, i+1]) / dx^2 + mu * (u_old[j+1, i] - 2*u_old[j, i] + u_old[j-1, i]) / dy^2

            # Equation 2.2 from chapter 5 of lecture notes
            u[j, i] = (rho * u_old[j, i] + dt * (diffusion + convection)) / rho
        end
    end
    
    # Solution for v-velocity 
    for i in 2:nx+1 # Iteration in x direction
        for j in 2:ny+1 # Iteration in y direction
            ve = 0.5 * (v_old[j, i] + v_old[j, i+1])
            vw = 0.5 * (v_old[j, i] + v_old[j, i-1])
            vn = 0.5 * (v_old[j, i] + v_old[j+1, i])
            vs = 0.5 * (v_old[j, i] + v_old[j-1, i])
            ue = 0.5 * (u_old[j-1, i+1] + u_old[j, i+1])
            uw = 0.5 * (u_old[j-1, i] + u_old[j, i])

            # Convective term 
            convection = -(rho * ve * ue - rho * vw * uw) / dx - (rho * vn^2 - rho * vs^2) / dy

            # Diffusive term 
            diffusion = mu * (v_old[j, i-1] - 2*v_old[j, i] + v_old[j, i+1]) / dx^2 + mu * (v_old[j+1, i] - 2*v_old[j, i] + v_old[j-1, i]) / dy^2

            # Equation 2.2 from chapter 5 of lecture notes
            v[j, i] = (rho * v_old[j, i] + dt * (diffusion + convection)) / rho
        end
    end
    return u, v 
end

function intermediate_velocity_upwind!(u, v, u_old, v_old, rho, mu, nx, ny, dx, dy, dt)
    #=
    Function to solve the conservation equation for the intermediate step (see lecture notes)

    Inputs:
    - u: u-velocity centered on the face - unknown 
    - v: v-velocity centered on the face - unknown
    - u_old: Known u-velocity centered on the face 
    - v_old: Known v-velocity centered on the face
    - rho: Fluid density [SI]
    - mu: Dynamic viscosity of fluid [SI]
    - nx: Number of internal elements in x 
    - ny: Number of internal elements in y 
    - dx: Step in x direction
    - dy: Step in y direction 
    - dt: Time step 

    Outputs:
    - u: u-velocity at intermediate step - Solution
    - v: v-velocity at intermediate step - Solution
    =#

    # Solution for u-velocity 
    for i in 3:nx+1 # Iteration in x direction
        for j in 2:ny+1 # Iteration in y direction
            ue = maximum([u_old[j,i], 0]) + minimum([0, u_old[j, i+1]])  #0.5 * (u_old[j, i]   + u_old[j, i+1])
            uw = maximum([u_old[j,i-1], 0]) + minimum([0, u_old[j, i]])   #0.5 * (u_old[j, i-1] + u_old[j, i])
            if ((v_old[j,i]) > 0.0 ) #if (v_old[j,i] > 0.0)
                un = u_old[j,i]
                us = u_old[j-1,i]
                vn = v_old[j,i]
                vs = v_old[j-1,i]
            else
                un = u_old[j+1, i]
                us = u_old[j,i]
                vn = v_old[j+1,i]
                vs = v_old[j,i]
            end

            # Convective term 
            convection = -(rho * ue^2 - rho * uw^2) / dx - (rho * vn * un - rho * vs * us) / dy

            # Diffusive term 
            diffusion = mu * (u_old[j, i-1] - 2*u_old[j, i] + u_old[j, i+1]) / dx^2 + mu * (u_old[j+1, i] - 2*u_old[j, i] + u_old[j-1, i]) / dy^2

            # Equation 2.2 from chapter 5 of lecture notes
            u[j, i] = (rho * u_old[j, i] + dt * (diffusion + convection)) / rho
        end
    end
    
    # Solution for v-velocity 
    for i in 2:nx+1 # Iteration in x direction
        for j in 3:ny+1 # Iteration in y direction
            vn = maximum([0, v_old[j,i]]) + minimum([0, v_old[j+1,i]])
            vs = maximum([0, v_old[j-1,i]]) + minimum([0, v_old[j,i]])

            if ((u_old[j,i]) > 0.0)
                vw = v_old[j-1,i]
                ve = v_old[j,i]
                ue = u_old[j,i]
                uw = u_old[j,i-1]
            else
                vw = v_old[j,i]
                ve = v_old[j,i+1]
                ue = u_old[j,i+1]
                uw = u_old[j,i]
            end

            # Convective term 
            convection = -(rho * ve * ue - rho * vw * uw) / dx - (rho * vn^2 - rho * vs^2) / dy

            # Diffusive term 
            diffusion = mu * (v_old[j, i-1] - 2*v_old[j, i] + v_old[j, i+1]) / dx^2 + mu * (v_old[j+1, i] - 2*v_old[j, i] + v_old[j-1, i]) / dy^2

            # Equation 2.2 from chapter 5 of lecture notes
            v[j, i] = (rho * v_old[j, i] + dt * (diffusion + convection)) / rho
        end
    end

    return u, v 
end


function rhie_chow_interpolation(u, v, p, rho, dx, dy, dt, nx, ny)
    # Initialize face velocity arrays
    u_face = zeros(Float64, ny+2, nx+3)  # For east/west faces
    v_face = zeros(Float64, ny+3, nx+2)  # For north/south faces
    
    # For east faces (u-velocity interpolation)
    for i in 3:nx+1
        for j in 2:ny+1
            # Standard linear interpolation
            avg_u = 0.5*(u[j,i] + u[j,i-1])
            
            # Pressure gradient term from momentum equation
            press_grad = dt/(2*rho*dx)*(p[j,i] - p[j,i-1])
            
            # Rhie-Chow correction term (third-order pressure dissipation)
            correction = dt/(2*rho*dx)*(p[j,i] - 2*p[j,i-1] + p[j,i-2])/4
            
            # Combined Rhie-Chow interpolated velocity
            u_face[j,i] = avg_u - press_grad + correction
        end
    end
    
    # For north faces (v-velocity interpolation)
    for i in 2:nx+1
        for j in 3:ny+1
            # Similar structure as u-velocity but in y-direction
            avg_v = 0.5*(v[j,i] + v[j-1,i])
            press_grad = dt/(2*rho*dy)*(p[j,i] - p[j-1,i])
            correction = dt/(2*rho*dy)*(p[j,i] - 2*p[j-1,i] + p[j-2,i])/4
            v_face[j,i] = avg_v - press_grad + correction
        end
    end
    
    return u_face, v_face
end

function accel(p, u, v, rho, nx, ny, dx, dy, dt)
  #=
  Function to calculate the velocity field u and v once the pressure is calculated with the
    intermediate velocity

    Inputs:
    - p: Pressure field centered in the cell 
    - u: u-velocity on cell face
    - v: v-velocity on cell face 
    - rho: Fluid density [SI]
    - dx: Step in x direction
    - dy: Step in y direction 
    - dt: Time step 

    Outputs:
    - u: Final corrected u-velocity field
    - v: Final corrected v-velocity field
  =#

    # Equation 2.3 from chapter 5 of lecture notes
    for i in 2:nx+1
        for j in 2:ny+1
            u[j, i] = (rho * u[j, i] - dt * (p[j, i] - p[j, i-1]) / (dx)) / rho
        end
    end

    # Equation 2.3 from chapter 5 of lecture notes
    for i in 2:nx+1
        for j in 2:ny+1
            v[j, i] = (rho * v[j, i] - dt * (p[j, i] - p[j-1, i]) / (dy)) / rho
        end
    end

    return u, v, u, v
end

function calc_divergence(u, v, nx, ny, dx, dy)
    #= 
    Function to calculate the divergence 

    Inputs:
    - u: u-velocity on cell face
    - v: v-velocity on cell face 
    - nx: Number of internal elements in x 
    - ny: Number of internal elements in y 
    - dx: Step in x direction
    - dy: Step in y direction 
    =#

    div = zeros(Float64, ny+2, nx+2)

    # Divergence calculation with forward difference 
    for i in 2:nx-1
        for j in 2:ny-1
            div[j, i] = (u[j, i+1] - u[j, i]) / dx + (v[j+1, i] - v[j, i]) / dy
        end
    end

    return div
end

function solve_pressure_poisson(u, v, p, nx, ny, dx, dy, dt, rho; TOL=1e-9, beta=1.872, maxIter = 300)
    #=
    Function to solve the Poisson equation using the 
    SOR (Successive Over-Relaxation) algorithm

    Inputs:
    - p: Pressure field centered in the cell 
    - u: u-velocity on cell face
    - v: v-velocity on cell face 
    - rho: Fluid density [SI]
    - dx: Step in x direction
    - dy: Step in y direction 
    - dt: Time step 
    - TOL: Accepted tolerance for iterative algorithm
    - beta: Relaxation coefficient 
    - maxIter: Maximum number of internal iterations for the iterative algorithm

    Outputs:
    - p: Pressure field solved by Poisson equation
    - iterations: Number of iterations needed to respect TOL tolerance
    =#

    # Get Rhie-Chow interpolated face velocities
    u_face, v_face = rhie_chow_interpolation(u, v, p, rho, dx, dy, dt, nx, ny)

    # Initializing pressure coefficients 
    a_e = ones(Float64, ny+2, nx+2) / (rho * dx^2)
    a_w = ones(Float64, ny+2, nx+2) / (rho * dx^2)
    a_n = ones(Float64, ny+2, nx+2) / (rho * dy^2)
    a_s = ones(Float64, ny+2, nx+2) / (rho * dy^2)
    
    # Setting pressure coefficients on boundaries to zero
    a_e[:, nx+1] .= 0.0
    a_w[:, 2]  .= 0.0
    a_n[ny+1, :] .= 0.0
    a_s[2, :]  .= 0.0
    
    # Central pressure coefficient
    a_p = -(a_e + a_w + a_n + a_s)
    
    maxError = 10.0
    iterations = 0
    
    while (abs(maxError) > TOL) && (iterations < maxIter)
        maxError = 0.0
        
        # Solving for pressure field 
        for i in 2:nx+1
            for j in 2:ny+1
                rhs = (u_face[j, i+1] - u_face[j, i]) / dx + (v_face[j+1, i] - v_face[j, i]) / dy
                rhs = rhs / dt - (a_w[j, i] * p[j, i-1] + a_e[j, i] * p[j, i+1] + 
                                 a_n[j, i] * p[j+1, i] + a_s[j, i] * p[j-1, i])
                
                p[j, i] = beta * rhs / a_p[j, i] + (1 - beta) * p[j, i]
            end
        end
        
        # Calculating residuals
        for i in 2:nx+1
            for j in 2:ny+1
                rhs = (u[j, i+1] - u[j, i]) / dx + (v[j+1, i] - v[j, i]) / dy
                error = (a_w[j, i] * p[j, i-1] + a_e[j, i] * p[j, i+1] + 
                         a_n[j, i] * p[j+1, i] + a_s[j, i] * p[j-1, i] + 
                         a_p[j, i] * p[j, i] - rhs / dt)
                
                if abs(error) > maxError
                    maxError = abs(error)
                end
            end
        end

        # Updating iteration count
        iterations += 1        
    end
    
    return p, iterations
end

function solve_cavity(domain, mesh, fluid_prop, bcs, solver)
    # Domain and mesh definition
    xe = domain.XE 
    ye = domain.YE 
    nx = mesh.NX 
    ny = mesh.NY
    dx, dy = xe / nx, ye / ny                       # mesh step-size

    x = collect(LinRange(0.5*dx, xe-0.5*dx, nx))    # x domain
    y = collect(LinRange(0.5*dy, ye-0.5*dy, ny))    # y domain 
    X = x' .* ones(length(y))                       # X - grid 
    Y = ones(length(x))'.* y                        # Y - grid 

    # Initialization 
    u, v, p, div, u_old, v_old = initialize_variables(nx, ny)

    # Input flow and problem parameters 
    rho = fluid_prop.RHO 
    mu = fluid_prop.MU 

    # Boundary Conditions 
    ut, ub = bcs.ut, bcs.ub
    vl, vr = bcs.vl, bcs.vr
    ul, ur = bcs.ul, bcs.ur
    vt, vb = bcs.vt, bcs.vb

    # Solver parameters 
    t = 0.0 

    ttake = @elapsed while (t < solver.END_TIME)
        # BC's 
        u, v, u_old, v_old = boundary_conditions!(u_old, v_old, ut, ub, ul, ur, vt, vb, vl, vr, nx, ny)
        
        # Time step update
        dt = stability_condition(u, v, dx, dy, CFL = solver.CFL)

        # Cálculo da velocidade intermediaria (passo estrela)
        u, v = intermediate_velocity!(u, v, u_old, v_old, rho, mu, nx, ny, dx, dy, dt)

        # Resolvendo a equação de Poisson no passo n+1
        p, _ = solve_pressure_poisson(u, v, p, nx, ny, dx, dy, dt, rho, maxIter=solver.MAXSOR_ITER, TOL = solver.SOR_TOL, beta = solver.SOR_BETA)

        # Cálculo da velocidade no passo n+1
        u, v, u_old, v_old = accel(p, u, v, rho, nx, ny, dx, dy, dt)

        t = t + dt 

        if (abs(maximum(u)) > 1e5)
            println("Divergence on u field")
            break
        end
    end
    println("Time to solve = ", ttake, "s")

    # Divergence calculation 
    div = calc_divergence(u, v, nx, ny, dx, dy)

    return x, y, X, Y, u, v, p, div
end

end