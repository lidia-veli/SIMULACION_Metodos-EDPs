'''
U_xx - (1+x)U_t + U = 0
'''

import numpy as np
import matplotlib.pyplot as plt


# PARÁMETROS
a, b = 0, 5  # x(t)
c, d = 0, 5  # t

N = 30
M = 30

h = (b-a)/N  
k = (d-c)/M  


# MATRIZ (t_j, x_i)
w = np.zeros((M+1, N+1))

# CONDICIONES CONTORNO

for j in range(M+1):  # eje t
    y_i = c + j*k
    w[j][0] = 0
    w[j][N] = 0

for i in range(N+1):  # eje x
    x_i = a + i*h
    w[0][i] = np.exp(-(x_i-2.5)**2)


# GAUSS-SEIDEL
for _ in range(100):  # sistema lineal
    for j in range(1,M):
        for i in range(1,N):
            # ej 1: U_xx -(1+x)U_t = 0
            #L = 2*k + h**2 + i*(h**3)
            #w[j][i] = (k/L)*(w[j][i+1]+w[j][i-1]) + ((h**2+i*(h**3))/L)*(w[j-1][i])

            # ej 2: U_xx - (1+x)U_t + U = 0
            L = 2*k + (h**2)*(1+i*h) - k*(h**2)
            w[j][i] = (k/L)*(w[j][i+1]+w[j][i-1]) + ((h**2)*(1+i*h)/L)*(w[j-1][i])

            # ej 3: U_xx - (1+x)U_t + xU = 0
            #L = 2*k + (h**2)*(1+i*h) - i*(h**3)*k
            #w[j][i] = (k/L)*(w[j][i+1]+w[j][i-1]) + ((h**2)*(1+h*i)/L)*(w[j-1][i])




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
ax.set_title('Solución EDP parabólica U_xx - (1+x)U_t + U = 0')
# ej 1: U_xx -(1+x)U_t = 0
# ej 2: U_xx - (1+x)U_t + U = 0
# ej 3: U_xx - (1+x)U_t + xU = 0
plt.show()
