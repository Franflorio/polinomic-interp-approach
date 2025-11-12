function yexp = interp_eval_point(xnodes, a, xexp, basis)
% INTERP_EVAL_POINT  Evalua un polinomio interpolante en un punto.
%   yexp = INTERP_EVAL_POINT(xnodes, a, xexp) asume coeficientes de NEWTON
%   (a0..an) obtenidos con NEWTON_DD y evalua en xexp usando Horner anidado.
%
%   yexp = INTERP_EVAL_POINT(xnodes, a, xexp, basis) permite elegir la base:
%     - basis = 'newton'   : a = [a0..an] de Newton, requiere xnodes
%     - basis = 'monomial' : a = [a0..an] monomicos ascendentes, xnodes no se usa
%
%   ENTRADAS
%     xnodes : nodos usados para construir 'a' si basis='newton'. Si
%              basis='monomial', puede pasarse [] o cualquier vector.
%     a      : vector de coeficientes (columna o fila).
%     xexp   : punto de evaluacion (escalar). Tambien admite vector.
%     basis  : 'newton' (default) o 'monomial'.
%
%   SALIDA
%     yexp   : valor p(xexp). Si xexp es vector, devuelve vector columna.
%
%   EJEMPLOS
%   % Newton (diferencias divididas):
%   x = [0;1;2]; y = [1;2;5];
%   [anew,~] = newton_dd(x, y);
%   y1 = interp_eval_point(x, anew, 1.5, 'newton');
%
%   % Monomica ascendente (p(x)=1 + 0*x + 1*x^2):
%   a = [1; 0; 1];
%   y2 = interp_eval_point([], a, 1.5, 'monomial');
%
%   Compatibilidad: MATLAB y Octave. ASCII-only.

if nargin < 4
    basis = 'newton';
end

a = a(:);
xq = xexp(:);
n = length(a) - 1;

if strcmpi(basis, 'newton')
    xnodes = xnodes(:);
    if length(xnodes) ~= n + 1
        error('Para basis=''newton'', length(xnodes) debe ser length(a).');
    end
    % Horner anidado para forma de Newton
    y = a(n+1) * ones(size(xq));
    k = n;
    while k >= 1
        y = a(k) + (xq - xnodes(k)) .* y;
        k = k - 1;
    end
elseif strcmpi(basis, 'monomial')
    % Horner clasico ascendente: p(x)=a0 + a1 x + ... + an x^n
    y = a(n+1) * ones(size(xq));
    k = n;
    while k >= 1
        y = a(k) + xq .* y;
        k = k - 1;
    end
else
    error('basis debe ser ''newton'' o ''monomial''.');
end

yexp = y;

% Si el usuario paso un escalar, opcionalmente devolver escalar
if isscalar(xexp)
    yexp = yexp(1);
end

end

