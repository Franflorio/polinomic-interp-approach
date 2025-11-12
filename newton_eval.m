function p = newton_eval(xnodes, a, xq)
% NEWTON_EVAL  Evalua el polinomio de Newton con coeficientes a.
%   p = NEWTON_EVAL(xnodes, a, xq) devuelve p(xq) usando Horner anidado.
%
%   ENTRADAS:
%     xnodes : nodos usados para construir a (mismo orden)
%     a      : [a0 a1 ... an] tal como devuelve NEWTON_DD
%     xq     : escalar o vector de puntos a evaluar
%
%   SALIDA:
%     p      : vector columna con p(xq)
%
%   EJEMPLO:
%     x = [0;1;2]; y = [1;2;5];
%     [a,~] = newton_dd(x,y);
%     xq = linspace(0,2,200).';
%     p = newton_eval(x,a,xq);
%     plot(x,y,'o',xq,p,'-'), grid on
%
%   Compatibilidad: MATLAB y Octave. ASCII-only.

xnodes = xnodes(:);
xq = xq(:);
n = length(a) - 1;

p = a(n+1) * ones(size(xq));
k = n;
while k >= 1
    p = a(k) + (xq - xnodes(k)) .* p;
    k = k - 1;
end
end

