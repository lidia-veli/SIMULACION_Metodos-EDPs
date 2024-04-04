import numpy as np
import matplotlib.pyplot as plt


# PARÁMETROS
a, b = 0, 1  # x
c, d = 0, 1  # y
N = 30
M = 30


# Tamaño del paso
h = (b-a)/N  
k = (d-c)/M  


# Matriz soluciones
w = np.zeros((N+1, M+1))

# CONDICIONES CONTORNO

for i in range(N+1):
    x_i = a + i*h
    w[i][0] = np.exp(-(x_i-0.25)**2)

for j in range(M+1):
    y_i = c + j*k
    w[0][j] = 0
    w[N][j] = 0


# Gauss-Seidel
L = 2*(k**2+h**2)
for _ in range(100):  # sistema lineal
    for j in range(1,M):  # da igual i,j que j,i
        for i in range(1,N):
            w[i][j] = ( ((k**2)/L)*(w[i+1][j]+w[i-1][j]) + ((h*k)/(2*L))*(w[i+1][j+1]+w[i-1][j-1]-w[i-1][j+1]-w[i+1][j-1]) + ((h**2)/L)*(w[i][j+1]+w[i][j-1]) )


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
ax.set_title('Solución EDP u_xx + 2u_xy + u_yy= 0')

plt.show()
