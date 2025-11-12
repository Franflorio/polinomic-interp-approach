function p = lagrange_eval_basic(xnodes, ynodes, xq)
% LAGRANGE_EVAL_BASIC  Interpolacion polinomica por Lagrange (forma clasica).
%   p = LAGRANGE_EVAL_BASIC(xnodes, ynodes, xq) evalua el polinomio
%   interpolante de Lagrange que pasa por los puntos (xnodes(i), ynodes(i))
%   en los puntos de consulta xq. Implementa la forma clasica de Lagrange:
%       p(x) = sum_{i=0..n} y_i * prod_{j!=i} (x - x_j) / (x_i - x_j)
%
%   REQUISITOS:
%     - xnodes y ynodes deben tener la misma longitud n+1.
%     - Todos los valores de xnodes deben ser distintos (nodos unicos).
%
%   ENTRADAS:
%     xnodes : vector de nodos (longitud n+1), no necesariamente ordenado.
%     ynodes : vector de valores f(xnodes) (longitud n+1).
%     xq     : vector de puntos a evaluar (m valores).
%
%   SALIDAS:
%     p      : vector columna de longitud m con p(xq(k)).
%
%   NOTAS NUMERICAS:
%     - Complejidad O(m * n^2). Adecuado para n moderado.
%     - Si xq coincide exactamente con un nodo xnodes(i), se devuelve
%       exactamente ynodes(i) (sin redondeo).
%     - Para mayor estabilidad y costo O(m * n), se recomienda la forma
%       baricentrica. Ver: BARYCENTRIC_WEIGHTS y BARYCENTRIC_EVAL.
%
%   EJEMPLO:
%     xnodes = [0; 1; 2];
%     ynodes = [1; 2; 5];
%     xq = linspace(0, 2, 9).';
%     p = lagrange_eval_basic(xnodes, ynodes, xq);
%     figure; plot(xnodes, ynodes, 'o', xq, p, '-'); grid on
%     legend('datos','Interpolante de Lagrange','Location','Best')
%     title('Interpolacion polinomica (Lagrange clasico)')
%
%   VEASE TAMBIEN:
%     BARYCENTRIC_WEIGHTS, BARYCENTRIC_EVAL, VANDERMONDE_INTERP,
%     POLY_EVAL_ASC, NEWTON_EVAL, CHEBYSHEV_NODES.
%
%   Compatibilidad: MATLAB y Octave. ASCII-only. Sin toolboxes.

% Normalizar a vectores columna
xnodes = xnodes(:);
ynodes = ynodes(:);
xq = xq(:);

% Chequeos basicos
n = length(xnodes) - 1;
if length(ynodes) ~= n + 1
    error('xnodes y ynodes deben tener la misma longitud (n+1 puntos)');
end
if length(unique(xnodes)) ~= length(xnodes)
    error('Los nodos en xnodes deben ser todos distintos');
end

m = length(xq);
p = zeros(m, 1);

k = 1;
while k <= m
    x = xq(k);

    % Si x coincide con un nodo, devolver exactamente ese valor
    idx = find(x == xnodes, 1);
    if ~isempty(idx)
        p(k) = ynodes(idx);
    else
        % Forma clasica: construir cada l_i(x) y acumular
        s = 0;
        i = 0;
        while i <= n
            num = 1;
            den = 1;
            j = 0;
            while j <= n
                if j ~= i
                    num = num * (x - xnodes(j+1));
                    den = den * (xnodes(i+1) - xnodes(j+1));
                end
                j = j + 1;
            end
            s = s + ynodes(i+1) * (num / den);
            i = i + 1;
        end
        p(k) = s;
    end

    k = k + 1;
end

end

