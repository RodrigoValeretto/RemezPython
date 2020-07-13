import numpy as np
from numpy.fft import fft
import matplotlib.pyplot as plt


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


def defineX(m, k, w, s, W, D):
    # Transforma os vetores em matrizes
    matrixm = m.reshape(1, len(m))
    matrixw = w[k].reshape(len(k), 1)
    # Multiplica as matrizes
    matrixwm = matrixw*matrixm
    # Calcula a matriz de cos
    coswm = np.cos(matrixwm)
    Wnp = np.array(W)

    # Calcula matriz que deve ser concat
    matrixsW = s/Wnp[k]

    # Concatena matrizes
    mF = np.zeros((coswm.shape[0], coswm.shape[1]+1))
    mF[:, :-1] = coswm
    mF[:, -1:] = matrixsW.reshape(coswm.shape[0], 1)

    # Calcula solução linear para mF \ D[k]
    matrixD = np.array(D)
    x = np.linalg.solve(mF, matrixD[k])
    return x


def defineA(x, M):
    a = []
    for i in range(0, int(M)+1):
        a.append(x[i])
    return a


def defineH(a, M):
    vet = []
    for i in range(int(M), 0, -1):
        vet.append(a[i])
    vet.append(2*a[0])
    for i in range(1, int(M)+1):
        vet.append(a[i])
    vet = np.array(vet)
    vet = vet/2
    return vet


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


def defineNovoK(vet1, vet2):
    res = []
    for i in range(len(vet1)):
        res.append(vet1[i])
    for i in range(len(vet2)):
        res.append(vet2[i])
    res = sorted(res)

    return res


def defineErroK(A, D, W, novok):
    errok = []
    for i in range(len(novok)):
        errok.append((A[novok[i]] - D[novok[i]])*W[novok[i]])
    
    return errok


def defineV(errok, delta, stop):
    res = []
    for i in range(0, len(errok)):
        if abs(errok[i]) >= abs(delta-stop):
            res.append(True)
    return res

def firamp(h, tipo, L):
    N = len(h)
    H = fft(h, 2*L)
    H = H[0:L+1]
    w = np.arange(0, L+1)*np.pi/L
    M = (N-1)/2

    if tipo == 1 or tipo == 2:
        H = np.exp(M*1j*w)*H
    else:
        H = -1j*np.exp(M*1j*w)*H

    A = H.real
    return A


def localMax(vet):
    n = len(vet)
    vcomp1 = []
    vcomp2 = []
    res = []

    vcomp1.append(vet[0]-1)
    for i in range(0, n-1):
        vcomp1.append(vet[i])

    for i in range(1, n):
        vcomp2.append(vet[i])
    vcomp2.append(vet[n-1]-1)

    for i in range(0, n):
        if vet[i] > vcomp1[i] and vet[i] > vcomp2[i]:
            res.append(i)

    return res


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

    while 1:
        # Resolvendo problema de interpolação
        x = defineX(m, k, w, s, W, D)
        a = defineA(x, M)
        delta = x[int(M)+1]
        h = defineH(a, M)
        A = firamp(h, 1, L)
        erro = (A - D)*W  # Erro do peso

        # Atualizando conjunto de referências
        novok = defineNovoK(localMax(erro), localMax(-erro))
        errok = defineErroK(A, D, W, novok)

        # Remove frequencias em que o erro do peso é menor que delta
        v = defineV(errok, delta, stop)
        print(v)
        
        break

    return 0  # h, delta


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
