function [h,delta] = remez(N,D,W)
% h = remez(N,D,W)
% Design de filtros tipo 1
%
% h : Resposta ao impulso de tamanho N
% N : Tamanho do filtro (Impar)
% D : Resposta Ideal (grid uniforme)
% W : Função de peso (grid uniforme)
% O tamanho de N deve ser igual ao de D

W = W(:);
D = D(:);
L = length(W)-1;
stop = 1e-8; % Numero pequeno definido para critério de parada
M = (N-1)/2;
R = M + 2; % R = Tamanho do conjunto de referencia

% Inicializa conjunto de referência (Espaçado igualmente aproximadamente onde W>0)
f = find(W>stop);
k = f(round(linspace(1,length(f),R)));
w = (0:L)'*pi/L;
m = 0:M;
s = (-1).^(1:R)'; % Sinais
while 1
% Resolvendo problema de interpolação
x = [cos(w(k)*m), s./W(k)] \ D(k);
a = x(1:M+1); 
delta = x(M+2); 
h = [a(M+1:-1:2); 2*a(1); a(2:M+1)]/2;
A = firamp(h,1,L)';
erro = (A-D).*W; % Erro do peso

% Atualizando conjunto de referências
novok = sort([localMax(erro); localMax(-erro)]);
errok = (A(novok)-D(novok)).*W(novok);

% Remove frequencias em que o erro do peso é menor que delta
v = abs(errok) >= (abs(delta)-stop);
novok = novok(v);
errok = errok(v);

% Garante a propriedade de alternação
v = gpalt(errok);
novok = novok(v);
errok = errok(v);

% Se o novo conjunto for muito grande, remove pontos até atingir tamanho
% adequado
while length(novok) > R
if abs(errok(1)) < abs(errok(length(novok)))
novok(1) = [];
else
novok(length(novok)) = [];
end
end

% Verifica se convergiu
if (max(errok)-abs(delta))/abs(delta) < stop
disp('Convergiu')
break
end
k = novok;
end
delta = abs(delta);
h = [a(M+1:-1:2); 2*a(1); a(2:M+1)]/2;
