'''
EJEMPLO: ecuación de ondas (DIF. PROGRESIVAS)
U_xx = (1/v^2)U_tt

condiciones:
* u(t,0) = u(t,L) = 0
* u(0,x) = f(x)
* U_t(0,x) = g(x)

f(x) = e^(-(x-2.5)**2)
g(x) = 0
'''

import numpy as np
import matplotlib.pyplot as plt
import sys

def verificar_velocidad_maxima(v, v_max):
    if v > v_max:
        print(f'v={v} debe ser menor que v_max={v_max}')
        sys.exit(1)
    return v



# PARÁMETROS
a, b = 0, 5  # x(t)
c, d = 0, 10  # t
N = 40
M = 400

# paso
h = (b-a)/N  # x(t)
k = (d-c)/M  # t

v = 0.5
# Condición de estabilidad (velocidad maxima)
v_max = h/k
v = verificar_velocidad_maxima(v, v_max)

# MATRIZ SOLUCIONES (t_j, x_i)
w = np.zeros((M+1, N+1))


# CONDICIONES FRONTERA

for j in range(1, M):  # eje t
    y_i = c + j*k
    w[j][0] = 0  # u(t,0)=0
    w[j][N] = 0  # u(t,L)=0


def f(x):
    return x*(b-x)
def g(x):
    return 0

for i in range(1, N):  # eje x
    x_i = a + i*h
    w[0][i] = f(x_i)  # u(0,x)=f(x)
    w[1][i] = w[0][i] + k*g(x_i)  # u(k,x) = u(0,x) + k*u_t(0,x)


# Método de EVOLUCIÓN (progresiva)
L = (k**2)*(v**2)/(h**2)
for j in range(1, M):  # DESDE 1
    for i in range(1, N):
        w[j+1][i]= (2-2*L)*w[j][i] + L*(w[j][i+1]+w[j][i-1]) - w[j-1][i]


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
ax.set_title('Solución Ec. Ondas (Met. Evolución)')

plt.show()
