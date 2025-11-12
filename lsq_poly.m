function [a, stats] = lsq_poly(x, y, n, method)
% LSQ_POLY  Ajuste polinómico por mínimos cuadrados (estándar).
%   [a,stats] = LSQ_POLY(x,y,n) devuelve los coeficientes 'a' (a0..an) del
%   polinomio de grado n que minimiza sum_i (y_i - p(x_i))^2 usando QR/backslash.
%   method opcional: 'qr' (default) o 'normal' (ecuaciones normales).
%
%   Referencia principal (internet):
%   Cleve Moler, Numerical Computing with MATLAB, cap. 5 (Least Squares).
%   Allí se describe el modelo y = X*beta y la solución con beta = X\y.
%
%   Ejemplo rápido:
%     x = [0; 2; 3; 5];  y = [-1; 0; 2; 1];  n = 2;
%     [a,stats] = lsq_poly(x,y,n);         % a = [a0;a1;a2]
%     yhat = polyval(flipud(a), x);        % flipud porque a es ascendente
%
%   Compatible con MATLAB y Octave. Sin toolboxes.

if nargin < 4
    method = 'qr';
end

x = x(:);
y = y(:);
m = length(x);

if length(y) ~= m
    error('x e y deben tener la misma longitud');
end
if n < 0
    error('El grado n debe ser >= 0');
end

% Matriz de diseño (Vandermonde ascendente): [1, x, x.^2, ..., x.^n]
X = ones(m, n+1);
j = 1;
while j <= n
    X(:, j+1) = X(:, j) .* x;
    j = j + 1;
end

% Resolver el problema LS
if strcmp(method,'qr')
    % Backslash: solución de mínimos cuadrados estable (QR/SVD bajo el capó)
    a = X \ y;
elseif strcmp(method,'normal')
    % Ecuaciones normales (menos estable si X es mal condicionada)
    G = X.' * X;
    c = X.' * y;
    a = G \ c;
else
    error('method debe ser ''qr'' o ''normal''');
end

% Métricas de ajuste
r = y - X*a;                 % residuos
SSE = sum(r.^2);             % suma de cuadrados de residuos
p = n + 1;
if m > p
    varianza = SSE / (m - p);
else
    varianza = NaN;
end
RMS = sqrt(SSE / m);

stats.SSE = SSE;
stats.RMS = RMS;
stats.varianza = varianza;
stats.residuos = r;

end

