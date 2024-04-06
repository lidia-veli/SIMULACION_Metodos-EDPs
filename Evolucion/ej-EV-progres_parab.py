'''
---------------------------------
---------------------------------
EJEMPLO: ecuación de calor (DIF. PROGRESIVAS)
v^2*U_xx = U_t

condiciones:
* u(t,0) = u(t,L) = 0
* u(0,x) = e^(-(x-2.5)**2)

'''
import numpy as np
import matplotlib.pyplot as plt


# PARÁMETROS
a, b = 0, 5  # x(t)
c, d = 0, 10  # t
N = 40
M = 400

v = 0.5

# paso
h = (b-a)/N  # x(t)
k = (d-c)/M  # t


# MATRIZ SOLUCIONES (t_j, x_i)
w = np.zeros((M+1, N+1))


# CONDICIONES FRONTERA

for j in range(1, M):  # eje t
    y_i = c + j*k
    w[j][0] = 0  # u(t,0)=0
    w[j][N] = 0  # u(t,L)=0

for i in range(1, N):  # eje x
    x_i = a + i*h
    w[0][i] = np.exp(-(x_i-2.5)**2)  # u(0,x)=f(x)


# Método de EVOLUCIÓN (progresiva)
L = k*(v**2)/(h**2)
for j in range(0, M):  # DESDE 0
    for i in range(1, N):
        w[j+1][i]= (1-2*L)*w[j][i] + L*(w[j][i+1]+w[j][i-1])


# Malla -----------------------------------------------------------------------
x = np.linspace(a, b, N+1)
t = np.linspace(c, d, M+1)
X, T = np.meshgrid(x, t)

# Creación de la figura 3D y los ejes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Graficación de la superficie
ax.plot_surface(X, T, w, cmap='viridis', edgecolor='none')
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('U(t.x)')
ax.set_title('Solución Ec. Calor (Met. Evolución)')

plt.show()
