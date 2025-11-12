function x = chebyshev_nodes(a, b, n)
% CHEBYSHEV_NODES  Nodos de Chebyshev de primer tipo en [a,b].
%   x = CHEBYSHEV_NODES(a, b, n) devuelve n+1 nodos recomendados para
%   interpolacion global estable (opcional para tu TP).
%
%   ENTRADAS
%     a,b : extremos del intervalo
%     n   : grado (entrega n+1 nodos)
%
%   SALIDA
%     x   : vector columna con n+1 nodos en [a,b]
%
%   Nota: utilidad estandar; no figura explicitamente en tus PDFs.
%   Compatibilidad: MATLAB y Octave. ASCII-only.

k = (0:n).';
x = (a + b)/2 + (b - a)/2 * cos( (2*k + 1) * pi / (2*(n+1)) );
end
