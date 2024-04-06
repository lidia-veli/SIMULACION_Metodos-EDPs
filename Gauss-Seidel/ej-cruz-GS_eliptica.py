'''
U_xx + U_xy + U_yy = 0
0<x<1, 0<y<1

condiciones:
u(0,y) = u(1,y) = u(x,1) = 0
u(x,0) = e^(-(x-0.25)^2)
'''

import numpy as np
import matplotlib.pyplot as plt


# PARÁMETROS
a, b = 0, 1  # x
c, d = 0, 1  # y(x)
N = 30
M = 30

# paso
h = (b-a)/N  
k = (d-c)/M  


# Matriz soluciones (x_i, y_j)
w = np.zeros((N+1, M+1))

# CONDICIONES CONTORNO

for i in range(N+1):  # eje x
    x_i = a + i*h
    w[i][0] = np.exp(-(x_i-0.25)**2)
    w[i][M] = 0

for j in range(M+1):  # eje y
    y_j = c + j*k
    w[0][j] = 0
    w[N][j] = 0


# GAUSS-SEIDEL
L = 2*(h**2 + k**2)
for _ in range(100):
    for j in range(1,M):  # da igual i,j que j,i
        for i in range(1,N):
            w[i][j] = ((k**2)/L)*(w[i+1][j]+w[i-1][j]) + ((h**2)/L)*(w[i][j+1]+w[i][j-1]) + ((h*k)/(4*L))*(w[i+1][j+1]+w[i-1][j-1]-w[i-1][j+1]-w[i+1][j-1])

# Malla
x = np.linspace(a, b, N+1)
y = np.linspace(c, d, M+1)
X, Y = np.meshgrid(x, y)

# Gráfica
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, w, cmap='viridis', edgecolor='none')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('U(x,y)')
ax.set_title('Solución EDP u_xx + u_xy + u_yy= 0')

plt.show()
