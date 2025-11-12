function [a, dd] = newton_dd(x, y)
% NEWTON_DD  Coeficientes de Newton via diferencias divididas.
%   [a, dd] = NEWTON_DD(x, y) calcula:
%     - a : vector fila [a0 a1 ... an] con a_k = f[x0,...,xk]
%     - dd: tabla de diferencias divididas (n x n), con dd(:,1)=y
%
%   ENTRADAS:
%     x, y  vectores (fila o columna) con n puntos (n >= 1), x distintos.
%
%   USO TIPICO:
%     [a,dd] = newton_dd(x,y);
%     % evaluar en xq con newton_eval(x,a,xq)
%
%   NOTAS:
%   - El orden de x fija la base (x-x0),(x-x0)(x-x1),...; mantenerlo
%     cuando evalues con NEWTON_EVAL.
%   - Si x tiene nodos muy cercanos o n grande, puede haber errores
%     numericos como en cualquier interpolacion global.
%
%   Compatibilidad: MATLAB y Octave. ASCII-only.

x = x(:);
y = y(:);

n = length(x);
if length(y) ~= n
    error('x e y deben tener la misma longitud');
end
if length(unique(x)) ~= n
    error('Los valores de x deben ser todos distintos');
end

dd = zeros(n, n);
dd(:,1) = y;

j = 2;
while j <= n
    i = 1;
    while i <= n - j + 1
        dd(i,j) = (dd(i+1,j-1) - dd(i,j-1)) / (x(i+j-1) - x(i));
        i = i + 1;
    end
    j = j + 1;
end

a = dd(1,1:n);
end

