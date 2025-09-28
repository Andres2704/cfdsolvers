using DelimitedFiles

# Including cavity_functions module 
include("./twoD_enthalpy.jl"); using .twoD_enthalpy


# Domain definition 
domain = (  
    L  = 1.0,                   # L length
    hi = 1.0,                    # Ice/water height 
    ha = 0.2                    # Solid part height
)

# Structured mesh properties
mesh = (
    NX = 240,                    # Number of x divisions for structured mesh 
    NZ = 240                     # Number of y divisions for structured mesh 
)

# Material properties 
mat_prop = (
    rhoa = 1000,                    # a-domain density [kg/m³]
    cp_a = 4185,                     # a-domain heat Capacity [J/(kg.K)]
    cp_i = 2060,                     # ice heat Capacity [J/(kg.K)]
    cp_w = 4185,                     # water heat Capacity [J/(kg.K)]
    rhoe = 1000,                    # e-domain density [kg/m³]
    l_a  = 0.6,                    # a-domain thermal conductivity [W/(m.K)]
    l_i  = 2.1,                    # Ice thermal conductivity [W/(m.K)]
    l_w  = 0.6,                    # Water thermal conductivity [W/(m.K)]
    T_m  = 273.15,                  # Melting temperature of water [K]
    Lm   = 333550,                  # Latent heat of water [J/kg]
)

# Boundary conditions 
bcs = ( 
    Trec = 214.12, 
    htc = 0
)

ics = (
    Ta_0 = 283.15,
    Te_0 = 263.15,
    S = 0.0 
)

# Solver properties 
solver = (
    CFL = 0.01,                  # CFL condition 
    END_TIME = 10,             # End time of simulation [s]
    MAXSOR_ITER = 300,          # Max SOR algorithm iterations 
    SOR_TOL = 1e-9,             # SOR algorithm tolerance 
    SOR_BETA = 1.872            # SOR beta parameter 
)

x, z, X, Z, H, T, rho, S, lambda, front, phi = solve_enthalpy(domain, mesh, mat_prop, bcs, ics, solver)


folder = "results"

# Check and create the folder if needed
if !isdir(folder)
    mkpath(folder)  # mkpath creates all necessary parent directories
end

writedlm(folder*"/enthalpy.txt", H, ',')
writedlm(folder*"/temperature.txt", T, ',')
writedlm(folder*"/density.txt", rho, ',')
writedlm(folder*"/source_term.txt", S, ',')
writedlm(folder*"/lambda.txt", lambda, ',')
writedlm(folder*"/X_coord.txt", X, ',')
writedlm(folder*"/Z_coord.txt", Z, ',')
writedlm(folder*"/front.txt", front, ',')
writedlm(folder*"/phi.txt", phi, ',')