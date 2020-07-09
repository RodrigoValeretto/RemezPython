function maxind = localMax(x)
% function maxind = localMax(x)
% Acha o indice do maximo local do vetor x
x = x(:);
n = length(x);
maxind = find(x > [x(1)-1;x(1:n-1)] & x > [x(2:n);x(n)-1]);