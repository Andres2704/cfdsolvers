using DelimitedFiles

# Including cavity_functions module 
include("./cavity_functions_rhiechow.jl"); using .cavity_functions_rhiechow


# Domain definition 
domain = (  
    XE = 1.0,                   # x-dimension in m
    YE = 1.0                    # y-dimension in m 
)

# Structured mesh properties
mesh = (
    NX = 16,                    # Number of x divisions for structured mesh 
    NY = 16                     # Number of y divisions for structured mesh 
)

# Fluid properties 
fluid_prop = (
    RHO = 1.0,                  # Fluid density [kg.m^-3]
    MU  = 0.01                  # Dynamic viscosity [Pa.s]
)

# Boundary conditions 
bcs = ( 
    ut = 1.0, ub = 0.0,         # Top and bottom u-velocity 
    vl = 0.0, vr = 0.0,         # Left and right v-velocity
    ul = 0.0, ur = 0.0,         # Left and right u-velocity
    vt = 0.0, vb = 0.0          # Top and bottom v-velocity
)

# Solver properties 
solver = (
    CFL = 0.1,                  # CFL condition 
    END_TIME = 4.0,             # End time of simulation [s]
    MAXSOR_ITER = 300,          # Max SOR algorithm iterations 
    SOR_TOL = 1e-9,             # SOR algorithm tolerance 
    SOR_BETA = 1.872            # SOR beta parameter 
)

x, y, X, Y, u, v, p, div = solve_cavity(domain, mesh, fluid_prop, bcs, solver)

folder = "results"

# Check and create the folder if needed
if !isdir(folder)
    mkpath(folder)  # mkpath creates all necessary parent directories
end

writedlm(folder*"/u_velocity_field_rhiechow.txt", u, ',')
writedlm(folder*"/v_velocity_field_rhiechow.txt", v, ',')
writedlm(folder*"/pressure_field_rhiechow.txt", p, ',')
writedlm(folder*"/divergence_field_rhiechow.txt", div, ',')
writedlm(folder*"/X_coordinate_field_rhiechow.txt", X, ',')
writedlm(folder*"/y_coordinate_field_rhiechow.txt", Y, ',')
writedlm(folder*"/x_coordinate_rhiechow.txt", x, ',')
writedlm(folder*"/y_coordinate_rhiechow.txt", y, ',')