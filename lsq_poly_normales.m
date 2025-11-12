function [a, stats] = lsq_poly_normales(x, y, n)
% LSQ_POLY_NORMALES  Ajuste por minimos cuadrados via ecuaciones normales.
%   [a,stats] = LSQ_POLY_NORMALES(x,y,n) calcula el polinomio
%   p(x) = a0 + a1 x + ... + an x^n que minimiza sum_k (y_k - p(x_k))^2
%   construyendo el sistema normal:
%       sum_{j=0..n} a_j * (sum_k x_k^(i+j)) = sum_k y_k * x_k^i,  i=0..n
%   (ver derivacion y forma matricial en los apuntes de catedra).
%
%   ENTRADAS
%     x, y : vectores de datos (misma longitud M)
%     n    : grado del polinomio (n >= 0)
%
%   SALIDAS
%     a      : coeficientes ascendentes [a0; a1; ...; an]
%     stats  : estructura con metricas
%              .M       = cantidad de puntos
%              .grado   = n
%              .SSE     = sum of squared errors
%              .RMS     = sqrt(SSE/M)   (valor medio cuadratico)
%              .sigma2  = SSE / (M - n - 1)   (varianza de la regresion)
%              .residuos= vector r = y - p(x)
%
%   Notas:
%   - Implementa textualmente las ecuaciones normales por sumas de potencias.
%   - La formula de RMS y la varianza coinciden con las de los apuntes.
%
%   Compatible con MATLAB y Octave. ASCII-only. Sin toolboxes.

x = x(:);
y = y(:);
M = length(x);

if length(y) ~= M
    error('x e y deben tener la misma longitud');
end
if n < 0
    error('El grado n debe ser >= 0');
end
if M == 0
    error('No hay datos');
end

% Precalcular sumas de potencias S(p) = sum_k x_k^p para p = 0..2n
S = zeros(2*n + 1, 1);
p = 0;
while p <= 2*n
    S(p+1) = sum(x.^p);
    p = p + 1;
end

% Vector b(i+1) = sum_k y_k * x_k^i, i = 0..n
b = zeros(n+1, 1);
i = 0;
while i <= n
    b(i+1) = sum(y .* (x.^i));
    i = i + 1;
end

% Matriz G(i+1,j+1) = sum_k x_k^(i+j) = S(i+j+1)
G = zeros(n+1, n+1);
i = 0;
while i <= n
    j = 0;
    while j <= n
        G(i+1, j+1) = S(i + j + 1);
        j = j + 1;
    end
    i = i + 1;
end

% Resolver ecuaciones normales
a = G \ b;

% Evaluar p(x) para metricas: construir matriz de diseno Vandermonde ascendente
X = ones(M, n+1);
j = 1;
while j <= n
    X(:, j+1) = X(:, j) .* x;
    j = j + 1;
end
yhat = X * a;
r = y - yhat;
SSE = sum(r.^2);

% Metricas (RMS y varianza como en la catedra)
RMS = sqrt(SSE / M);
den = M - n - 1;
if den > 0
    sigma2 = SSE / den;
else
    sigma2 = NaN;
end

stats.M = M;
stats.grado = n;
stats.SSE = SSE;
stats.RMS = RMS;
stats.sigma2 = sigma2;
stats.residuos = r;

end

