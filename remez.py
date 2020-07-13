import numpy as np
from numpy.fft import fft
import matplotlib.pyplot as plt


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

    len(s)
    len(Wnp)
    len(Wnp[k])
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


def defineV(errok, delta, stop):
    res = []
    for i in range(0, len(errok)):
        if abs(errok[i]) >= (abs(delta)-stop):
            res.append(True)
        else:
            res.append(False)
    return res


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


def atualizaNEK(nek, v, inteiro):
    res = []
    if not inteiro:
        for i in range(len(v)):
            if v[i]:
                res.append(nek[i])
    else:
        for i in range(len(v)):
            res.append(nek[v[i]])

    return res


def ajustaTamConj(novok, errok, R):
    while len(novok) > R:
        if abs(errok[0]) < abs(errok[len(novok)-1]):
            novok.remove(novok[0])
        else:
            novok.remove(novok[len(novok)-1])
    return novok


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


def gpalt(vet):
    res = []
    xe = vet[0]
    xv = 0
    for i in range(1, len(vet)):
        if np.sign(vet[i]) == np.sign(xe):
            if abs(vet[i]) > abs(xe):
                xe = vet[i]
                xv = i
        else:
            res.append(xv)
            xe = vet[i]
            xv = i

    res.append(xv)
    return res


def remez(N, D, W, itmax):
    # Função que implementa o algoritmo de remez
    # para construção de filtros FIR
    L = len(W) - 1
    stop = 10**(-8)
    M = int((N-1)/2)
    R = M + 2

    # Inicializando conjunto de referências
    f = find(W, stop, True)
    k = defineK(R, f)
    w = np.arange(0, L+1)*np.pi/L
    m = np.arange(0, M+1)
    s = defineS(R)

    it = 1
    while it <= itmax:
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
        novok = atualizaNEK(novok, v, False)
        errok = atualizaNEK(errok, v, False)

        # Garante propriedade de alternação
        v = gpalt(errok)
        novok = atualizaNEK(novok, v, True)
        errok = atualizaNEK(errok, v, True)

        # Se o novo conjunto for muito grande, remove pontos até atingir tamanho
        # adequado
        novok = ajustaTamConj(novok, errok, R)

        # Confere se o algoritmo convergiu
        if (max(errok)-abs(delta))/abs(delta) < stop:
            print("Convergiu com " + str(it) + " iterações.")
            break

        if it == itmax:
            print("Algoritmo não convergiu")

        k = novok
        it += 1

    delta = abs(delta)
    h = defineH(a, M)

    return h, A, delta



'''# Script de teste1
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

itmax = 1000  # numero máximo de iterações

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
h, A, delta = remez(N, D, W, itmax)
'''

'''# Script de teste2
itmax = 1000  # numero máximo de iterações
N = 31
Kp = 1
Ks = 4
wp = np.dot(0.26, np.pi)
ws = np.dot(0.34, np.pi)
wo = np.dot(0.3, np.pi)
L = 1000
w = np.dot(np.concatenate([np.arange(0, L+1)]), np.pi) / L
W = np.dot(Kp, (w <= wp)) + np.dot(Ks, (w >= ws))
D = 1*(w <= wo)
h, A, delta = remez(N, D, W, itmax)
'''

# Script de teste3
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


itmax = 1000  # numero máximo de iterações

N = 105
ws1 = 0.20*np.pi
wp1 = 0.25*np.pi
wp2 = 0.63*np.pi
ws2 = 0.68*np.pi
Ks1 = 10
Kp = 1
Ks2 = 1

wo = (ws1+wp1)/2
w1 = (wp1+ws2)/2

L = 4000
w = np.arange(0, L+1)*np.pi/L
W = defineW(w)
D = defineD(w)
h, A, delta = remez(N, D, W, itmax)


# Plot
n = np.arange(0, N)
fig, ax = plt.subplots(2)
ax[0].stem(h, use_line_collection=True)
ax[0].set_title("h[n] do Filtro")
ax[1].plot(A)
ax[1].set_title("Resposta em frequencia")
plt.show()
