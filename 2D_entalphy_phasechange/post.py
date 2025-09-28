import numpy as np 
import matplotlib.pyplot as plt 

folder = "results"
H = np.loadtxt(folder+"/enthalpy.txt", delimiter = ",")
rho = np.loadtxt(folder+"/density.txt", delimiter = ",")
lamb = np.loadtxt(folder+"/lambda.txt", delimiter = ",")
S = np.loadtxt(folder+"/source_term.txt", delimiter = ",")
T = np.loadtxt(folder+"/temperature.txt", delimiter = ",")
X = np.loadtxt(folder+"/X_coord.txt", delimiter = ",")
Z = np.loadtxt(folder+"/Z_coord.txt", delimiter = ",")
phi = np.loadtxt(folder+"/phi.txt", delimiter = ",")
front = np.loadtxt(folder+"/front.txt", delimiter = ",")

t_sim = np.linspace(0, 10, len(front))

t = np.linspace(0, 10, 1000)
xsi = 0.29573601432884933
lw = 0.6
rhow = 1000
cpw = 4185 
x_st = 2*xsi * np.sqrt(lw*t/(rhow*cpw))

# List of fields and titles
fields = [H, rho, lamb, S, T, phi]
titles = ["Enthalpy (H)", "Density (ρ)", "Thermal Conductivity (λ)", "Source Term (S)", "Temperature (T)", "Phi"]

# Create subplots
# fig, axes = plt.subplots(2, 3, figsize=(18, 10))
# axes = axes.flatten()

# for i, (field, title) in enumerate(zip(fields, titles)):
#     contour = axes[i].pcolormesh(X, Z, field, shading='auto', antialiased=False, cmap='viridis')
#     axes[i].set_title(title)
#     axes[i].set_xlabel("X")
#     axes[i].set_ylabel("Z")
#     fig.colorbar(contour, ax=axes[i], orientation='vertical')

# # Hide unused subplot if odd number
# if len(fields) < len(axes):
#     axes[-1].axis('off')

# plt.tight_layout()


plt.figure()
plt.plot(t, x_st*1000, color = 'b', marker = 's', markevery = 10, label = 'Analytical Solution')
plt.plot(t_sim, front*1000, color = 'r', marker = 'o', markevery = 10, label = 'Simulation')
plt.grid()
plt.xlabel('t [s]')
plt.ylabel('Melting front [mm]')

plt.figure()
plt.plot(Z[:,100], H[:,100], color = 'b')
plt.grid()
plt.xlabel('z [m]')
plt.ylabel('Enthalpy [J/kg]')

plt.show()