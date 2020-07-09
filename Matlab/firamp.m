function [A,w] = firamp(h,type,L)
% [A,w] = firamp(h,type,L)
% Resposta de amplitude de um filtro FIR de fase linear
% A : Resposta de amplitude
% w : grid de frequencia [0:L]*pi/L
% h : resposta ao impulso
% type : [1,2,3,4]
% L : Densidade de frequência(Opcional, default = 2^10)

h = h(:)';		% transforma h em um vetor de linha
N = length(h);  % tamanho de h
if nargin < 3
	L = 2^10;	% tamanho do grid
end
H = fft(h,2*L);		% fft
H = H(1:L+1);		% seleciona[0,pi]
w = [0:L]*pi/L;		% grid de frequencia
M = (N-1)/2;		
if (type == 1)|(type == 2)
  H = exp(M*j*w).*H;	% Tipo I e II
else
  H = -j*exp(M*j*w).*H;	% Tipo III e IV
end
A = real(H);		% Descarta a parte imaginaria zero


