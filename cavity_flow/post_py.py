import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
from matplotlib.colors import Normalize

# Carrega dados Ghia et al.
initial_velocity = 1.0 

folder = "results"
A = pd.read_csv("Re_100.csv", header=None).values
B = pd.read_csv("re_100y.csv", header=None).values
u = np.loadtxt(folder+"/u_velocity_field_staggered.txt", delimiter = ",")
v = np.loadtxt(folder+"/v_velocity_field_staggered.txt", delimiter = ",")
p = np.loadtxt(folder+"/pressure_field_staggered.txt", delimiter = ",")
x = np.loadtxt(folder+"/x_coordinate_staggered.txt", delimiter = ",")
y = np.loadtxt(folder+"/y_coordinate_staggered.txt", delimiter = ",")
X = np.loadtxt(folder+"/X_coordinate_field_staggered.txt", delimiter = ",")
Y = np.loadtxt(folder+"/Y_coordinate_field_staggered.txt", delimiter = ",")
div = np.loadtxt(folder+"/divergence_field_staggered.txt", delimiter = ",")

u_rc = np.loadtxt(folder+"/u_velocity_field_rhiechow.txt", delimiter = ",")
v_rc = np.loadtxt(folder+"/v_velocity_field_rhiechow.txt", delimiter = ",")
p_rc = np.loadtxt(folder+"/pressure_field_rhiechow.txt", delimiter = ",")
x_rc = np.loadtxt(folder+"/x_coordinate_rhiechow.txt", delimiter = ",")
y_rc = np.loadtxt(folder+"/y_coordinate_rhiechow.txt", delimiter = ",")
X_rc = np.loadtxt(folder+"/X_coordinate_field_rhiechow.txt", delimiter = ",")
Y_rc = np.loadtxt(folder+"/Y_coordinate_field_rhiechow.txt", delimiter = ",")
div_rc = np.loadtxt(folder+"/divergence_field_rhiechow.txt", delimiter = ",")

nx = len(x)  
ny = len(y)  
nx_rc = len(x_rc) 
ny_rc  = len(y_rc) 

# Interpola velocidades centradas
U = 0.5 * (u[1:-1,1:-1] + u[1:-1,2:])
V = 0.5 * (v[1:-1,1:-1] + v[2:,1:-1])

norm = Normalize(vmin=0, vmax=1)

V_mag = np.sqrt(U**2 + V**2)
V_mag_rc = np.sqrt(u_rc**2 + v_rc**2)

# Interpola perfis para comparação
Uq = np.interp(y, y, U[:, nx//2])
Vq = np.interp(x, x, V[ny//2, :])

# np.savetxt('cavity_solver/Uq_solver.dat', np.c_[y, Uq], delimiter=',')
# np.savetxt('cavity_solver/Vq_solver.dat', np.c_[x, Vq], delimiter=',')


# Gráficos de comparação
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(A[:,0], A[:,1], label='Ghia et al.')
plt.plot(y, Uq / initial_velocity, 'o', label='Staggered Grid')
plt.plot(y_rc, u_rc[1:-1, nx_rc//2] / initial_velocity, 'x', label='Rhie and Chow Interpolation', color = 'r')
plt.title('U @ x=0.5')
plt.xlabel('y')
plt.ylabel('U/U0')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(B[:,0], B[:,1], label='Ghia et al.')
plt.plot(x, Vq / initial_velocity, 'o', label='Staggered Grid')
plt.plot(x_rc, v_rc[ny_rc//2, 1:-1] / initial_velocity, 'x', label='Rhie and Chow Interpolation', color = 'r')
plt.title('V @ y=0.5')
plt.xlabel('x')
plt.ylabel('V/U0')
plt.legend()
plt.tight_layout()


# Vetores de velocidade
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.streamplot(X, Y, U, V, density = 2, color = V_mag, norm=norm)
plt.title("Staggered Grid")
plt.xlabel("x")
plt.ylabel("y")
plt.xlim([0,x[-1]])
plt.ylim([0,y[-1]])
plt.colorbar(label='Vmag')
plt.axis("equal")
plt.tight_layout()

plt.subplot(1, 2, 2)
plt.streamplot(X_rc, Y_rc, u_rc[1:-1,1:-1], v_rc[1:-1,1:-1], density = 2, color = V_mag_rc[1:-1,1:-1], norm=norm)
plt.title("Rhie & Chow Interpolation")
plt.xlabel("x")
plt.ylabel("y")
plt.xlim([0,x[-1]])
plt.ylim([0,y[-1]])
plt.colorbar(label='Vmag')
plt.axis("equal")
plt.tight_layout()

# # Campo de pressão
# plt.figure(figsize=(6, 5))
# plt.contourf(X, Y, p[1:-1,1:-1], levels=50, cmap='viridis')
# plt.colorbar(label='Pressão')
# plt.title("Campo de Pressão")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.axis("equal")
# plt.tight_layout()
    
# # Plot divergence
# plt.figure(figsize=(6, 5))
# plt.contourf(X, Y, div[1:ny+1, 1:nx+1], levels=50, cmap='coolwarm')
# plt.colorbar(label='Divergente')
# plt.title("Campo do Divergente")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.axis("equal")
# plt.tight_layout()
plt.show()