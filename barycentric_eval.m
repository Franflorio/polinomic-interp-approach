function p = barycentric_eval(xnodes, ynodes, w, xq)
% BARYCENTRIC_EVAL  Evalua el interpolante baricentrico en puntos xq.
%   p = BARYCENTRIC_EVAL(xnodes, ynodes, w, xq) devuelve p(xq) usando
%   la formula baricentrica estable. Si xq coincide con un nodo, devuelve
%   exactamente el y correspondiente.
%
%   ENTRADAS
%     xnodes : nodos (n+1)
%     ynodes : valores f(xnodes)
%     w      : pesos baricentricos (BARYCENTRIC_WEIGHTS)
%     xq     : puntos de evaluacion (m)
%
%   SALIDA
%     p      : valores p(xq) (vector columna m)
%
%   Compatibilidad: MATLAB y Octave. ASCII-only.

xnodes = xnodes(:);
ynodes = ynodes(:);
w = w(:);
xq = xq(:);

n = length(xnodes);
if length(ynodes) ~= n || length(w) ~= n
    error('xnodes, ynodes y w deben tener la misma longitud.');
end

m = length(xq);
p = zeros(m,1);

k = 1;
while k <= m
    xd = xq(k) - xnodes;
    idx = find(xd == 0, 1);
    if ~isempty(idx)
        p(k) = ynodes(idx);
    else
        num = sum( (w .* ynodes) ./ xd );
        den = sum( w ./ xd );
        p(k) = num / den;
    end
    k = k + 1;
end
end
