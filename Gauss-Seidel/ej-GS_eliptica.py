'''
Ec. Poisson
U_xx + U_yy = f(x,y)
0<x<1, 0<y<1
f(x,y) = x*y*(1-x)*(1-y)

condiciones contorno:
* u(0,y) = u(x,0) = u(x,1) = u(1,y)
'''

import numpy as np
import matplotlib.pyplot as plt


# PARÁMETROS
a, b = 0, 1  # x
c, d = 0, 1  # y(x)

N = 30
M = 30

h = (b-a)/N  
k = (d-c)/M  


# MATRIZ (x_i, y_i)
w = np.zeros((N+1, M+1))

# CONDICIONES CONTORNO

for i in range(N+1):  # eje x
    x_i = a + i*h
    w[i][0] = 0
    w[i][N] = 0

for j in range(M+1):  # eje y
    y_j = c + j*k
    w[0][j] = 0
    w[N][j] = 0


# GAUSS-SEIDEL
def f(i,j):
    x_i = a + i*h
    y_j = c + j*k
    return x_i*y_j*(1-x_i)*(1-y_j)

for _ in range(100):  # sistema lineal
    for j in range(1,M):
        for i in range(1,N):
            L = 2*(h**2 + k**2)
            w[i][j] = ((k**2)/L)*(w[i+1][j]+w[i-1][j]) + ((h**2)/L)*(w[i][j+1]+w[i][j-1]) - ((h**2)*(k**2)/L)*f(i,j)



# Malla --------------------------------------------------------------------------------------------------
x = np.linspace(a, b, N+1)
t = np.linspace(c, d, M+1)
X, T = np.meshgrid(x, t)

# Gráfica
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, T, w, cmap='viridis', edgecolor='none')
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('U(t,x)')
ax.set_title('Solución EDP elíptica U_xx + U_yy = x*y*(1-x)*(1-y)')
plt.show()
