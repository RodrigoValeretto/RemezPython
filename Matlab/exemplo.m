% Exemplo
N = 31; Kp = 1; Ks = 4;
wp = 0.26*pi; ws = 0.34*pi; wo = 0.3*pi;
L = 1000;
w = [0:L]*pi/L;
W = Kp*(w<=wp) + Ks*(w>=ws);
D = (w<=wo);
h = remez(N,D,W);