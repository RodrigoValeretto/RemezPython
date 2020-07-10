import numpy as np


def find(vet, comp, maior):
    res = []
    if maior:
        for i in range(0, len(vet)):
            if vet[i] >= comp:
                res.append(i)
    else:
        for i in range(0, len(vet)):
            if vet[i] <= comp:
                res.append(i)

    return res


def defineK(R, f):
    vet = np.linspace(0, len(f)-1, int(R))
    aux = []
    for i in range(0, len(vet)):
        aux.append(f[round(vet[i])])
    return aux

def defineS(R):
    s = []
    for i in range(1, int(R+1)):
        s.append(((-1)**i))
    return s

def defineX(m, k, w, s, D):
    print(np.asmatrix(m))


# Função que implementa o algoritmo de remez
# para construção de filtros FIR
def remez(N, D, W):
    L = len(W) - 1
    stop = 10**(-8)
    M = (N-1)/2
    R = M + 2

    # Inicializando conjunto de referências
    f = find(W, stop, True)
    k = defineK(R, f)
    w = np.arange(0, L+1)*np.pi/L
    m = np.arange(0, M+1)
    s = defineS(R)
    '''
    while 1:
        x = [cos(w(k)*m), s/W[k]] \ D[k]
    '''
    return 0  # h, delta


def defineW(w):
    vet = []
    for i in range(0, len(w)):
        vet.append(Ks1*(w[i] <= ws1) + Kp*((w[i] >= wp1)
                                           and (w[i] <= wp2)) + Ks2*(w[i] >= ws2))
    return vet


def defineD(w):
    vet = []
    for i in range(0, len(w)):
        vet.append(int((wp1 <= w[i]) and (w[i] <= wp2)))
    return vet


# Script de teste
N = 103
ws1 = 0.20*np.pi
wp1 = 0.25*np.pi
wp2 = 0.60*np.pi
ws2 = 0.70*np.pi
Ks1 = 10
Kp = 1
Ks2 = 1
wo = (ws1+wp1)/2
w1 = (wp1+ws2)/2
L = 4000
w = np.arange(0, L+1)*np.pi/L
W = defineW(w)
D = defineD(w)
h = remez(N, D, W)
