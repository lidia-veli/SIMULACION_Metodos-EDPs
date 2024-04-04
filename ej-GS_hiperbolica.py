import numpy as np
import matplotlib.pyplot as plt


# PARÁMETROS
a, b = 0, 0.5  # x
c, d = 0, 0.5  # y
N = 30
M = 30


# Tamaño del paso
h = (b-a)/N  
k = (d-c)/M  


# Matriz soluciones
w = np.zeros((N+1, M+1))

# CONDICIONES CONTORNO
def f(x):
    return np.exp(-(x-0.25)**2)

for i in range(N+1):
    x_i = a + i*h
    w[i][0] = f(x_i)

for j in range(M+1):
    y_i = c + j*k
    w[0][j] = 0
    w[N][j] = 0


# Gauss-Seidel
for _ in range(100):  # sistema lineal
    for i in range(1,N):
        for j in range(1,M):
            w[i][j] = (1/2)*(w[i+1][j] + w[i-1][j]) + (h/(8*k))*(w[i+1][j+1] + w[i-1][j-1] - w[i-1][j+1] - w[i+1][j-1])

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
ax.set_title('Solución EDP u_xx + u_xy = 0')

plt.show()
