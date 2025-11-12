function a = vandermonde_interp(x, y)
% VANDERMONDE_INTERP  Interpolacion polinomica via coeficientes indeterminados.
% Devuelve coeficientes ascendentes a0..an tales que p(xi)=yi.
% CUIDADO: para n grande puede ser mal condicionado.

x = x(:);
y = y(:);
n = length(x) - 1;
if length(y) ~= n + 1
    error('x e y deben tener la misma longitud (n+1 puntos)');
end

V = ones(n+1, n+1);
j = 1;
while j <= n
    V(:, j+1) = V(:, j) .* x;
    j = j + 1;
end

a = V \ y;
a = a(:);
end

