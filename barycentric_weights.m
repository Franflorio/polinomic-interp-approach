function w = barycentric_weights(x)
% BARYCENTRIC_WEIGHTS  Pesos baricentricos para interpolacion de Lagrange.
%   w = BARYCENTRIC_WEIGHTS(x) devuelve los pesos w_i = 1/prod_{j!=i}(x_i-x_j)
%   para los nodos x(:). Complejidad O(n^2). Adecuado para n moderado.
%
%   ENTRADA
%     x : vector de nodos (fila o columna), todos distintos.
%
%   SALIDA
%     w : vector columna de pesos baricentricos.
%
%   Compatibilidad: MATLAB y Octave. ASCII-only.

x = x(:);
n = length(x);
if length(unique(x)) ~= n
    error('Los nodos x deben ser todos distintos.');
end

w = ones(n,1);
i = 1;
while i <= n
    j = 1;
    while j <= n
        if j ~= i
            w(i) = w(i) / (x(i) - x(j));
        end
        j = j + 1;
    end
    i = i + 1;
end
end
