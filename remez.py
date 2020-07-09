import numpy as np
from random import *

N = 13                            # Filter Length
Kp = 1                            # Peso de passagem
Ks = 2                            # Peso de parada
wp = 0.4*np.pi                    # borda de passagem
ws = 0.5*np.pi                    # borda de parada
wo = (wp+ws)/2                    # freq de corte
L = 1000                          # tamanho do grid
w = np.linspace(0, 1, L)*np.pi/L  # Frequencia
D = (w <= wo)                     # função desejada
M = (N-1)/2
R = M + 2


def W(w):
    if w <= wp:
        return Kp*w
    elif w >= ws:
        return Ks*w
    elif w >= ws and w <= wp:
        return Kp*w + Ks*w
    else:
        return 0


k = []
while len(k) != 8:
    aux = randint(0, L)
    if W(w[aux]) != 0:
        k.append(aux)

ref = w[k]/np.pi

# Iteration 1
m = range(0, int(M))
s = []
for i in range(1, int(R+1)):
    s.append((-1)**i)