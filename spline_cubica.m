function S = spline_cubica(x, y, tipo, fp0, fpn)
% SPLINE_CUBICA  Construye un spline cubico (natural o sujeto/clamped).
%   S = SPLINE_CUBICA(x,y) construye el spline cubico NATURAL que interpola
%   los datos (x(i), y(i)), con x estrictamente creciente.
%
%   S = SPLINE_CUBICA(x,y,'natural') idem anterior.
%
%   S = SPLINE_CUBICA(x,y,'clamped',fp0,fpn) construye el spline SUJETO
%   especificando las derivadas en los extremos: S'(x0)=fp0, S'(xn)=fpn.
%
%   Salida S: estructura con campos
%     .x  : nodos x(1..n)
%     .a  : a(i) = y(i)                      (i = 1..n-1)
%     .b  : coeficientes lineales por tramo  (i = 1..n-1)
%     .c  : coeficientes cuadraticos por tramo (i = 1..n-1); c(n) no se usa
%     .d  : coeficientes cubicos por tramo   (i = 1..n-1)
%     .tipo, .fp0, .fpn : metadatos del borde
%
%   Notas:
%   - Implementa el sistema tridiagonal clasico para c(0..n) y luego
%     b(i) = (a(i+1)-a(i))/h(i) - h(i)*(2*c(i)+c(i+1))/3
%     d(i) = (c(i+1)-c(i)) / (3*h(i))
%   - Para 'clamped' se usan las ecuaciones de borde:
%       2*h0*c0 + h0*c1 = 3*((a1-a0)/h0 - fp0)
%       h_{n-1}*c_{n-1} + 2*h_{n-1}*c_n = 3*(fpn - (an-a_{n-1})/h_{n-1})
%
%   Compatibilidad: MATLAB y Octave.

if nargin < 3 || isempty(tipo)
    tipo = 'natural';
end

x = x(:); y = y(:);
n = length(x);
if length(y) ~= n
    error('x e y deben tener la misma longitud.');
end
if n < 2
    error('Se requieren al menos dos puntos.');
end
if any(diff(x) <= 0)
    error('x debe ser estrictamente creciente.');
end

% Preparativos
h = diff(x);          % h(i) = x(i+1) - x(i)
a = y;                % a(i) = y(i)

% Construir sistema tridiagonal para c0..c_n
T = zeros(n, n);
rhs = zeros(n, 1);

if strcmpi(tipo, 'natural')
    % ecuaciones interiores
    i = 2;
    while i <= n-1
        T(i, i-1) = h(i-1);
        T(i, i)   = 2*(h(i-1) + h(i));
        T(i, i+1) = h(i);
        rhs(i) = 3*((a(i+1)-a(i))/h(i) - (a(i)-a(i-1))/h(i-1));
        i = i + 1;
    end
    % bordes naturales
    T(1,1) = 1;       rhs(1) = 0;
    T(n,n) = 1;       rhs(n) = 0;

elseif strcmpi(tipo, 'clamped')
    if nargin < 5
        error('Para ''clamped'' debe especificar fp0 y fpn.');
    end
    % ecuacion izquierda
    T(1,1) = 2*h(1);
    T(1,2) = h(1);
    rhs(1) = 3*((a(2)-a(1))/h(1) - fp0);

    % ecuaciones interiores
    i = 2;
    while i <= n-1
        T(i, i-1) = h(i-1);
        T(i, i)   = 2*(h(i-1) + h(i));
        T(i, i+1) = h(i);
        rhs(i) = 3*((a(i+1)-a(i))/h(i) - (a(i)-a(i-1))/h(i-1));
        i = i + 1;
    end

    % ecuacion derecha
    T(n, n-1) = h(n-1);
    T(n, n)   = 2*h(n-1);
    rhs(n) = 3*(fpn - (a(n)-a(n-1))/h(n-1));

else
    error('tipo debe ser ''natural'' o ''clamped''.');
end

% Resolver para c(1..n) correspondiendo a c_0..c_n
c_all = T \ rhs;

% Calcular b(i), c(i), d(i) por tramo i=1..n-1
c = c_all(1:n-1);
b = zeros(n-1,1);
d = zeros(n-1,1);

i = 1;
while i <= n-1
    b(i) = (a(i+1)-a(i))/h(i) - h(i)*(2*c_all(i) + c_all(i+1))/3;
    d(i) = (c_all(i+1) - c_all(i)) / (3*h(i));
    i = i + 1;
end

% Guardar en estructura
S.x = x;
S.a = a(1:n-1);
S.b = b;
S.c = c;
S.d = d;
S.tipo = lower(tipo);
if strcmpi(tipo,'clamped')
    S.fp0 = fp0;
    S.fpn = fpn;
else
    S.fp0 = [];
    S.fpn = [];
end

end

