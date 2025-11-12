function v = spline_eval(S, xq)
% SPLINE_EVAL  Evalua el spline S en puntos xq.
%   v = SPLINE_EVAL(S, xq) devuelve S(xq). S es la estructura creada por
%   SPLINE_CUBICA. Para cada xq se localiza el tramo adecuado.
%
%   Compatibilidad: MATLAB y Octave.

xq = xq(:);
x = S.x;
a = S.a; b = S.b; c = S.c; d = S.d;
n = length(x);

v = zeros(size(xq));

k = 1;
while k <= length(xq)
    t = xq(k);
    if t <= x(1)
        i = 1;
    elseif t >= x(n)
        i = n-1;
    else
        i = find(t >= x(1:n-1) & t < x(2:n), 1, 'first');
    end
    dx = t - x(i);
    v(k) = a(i) + b(i)*dx + c(i)*dx^2 + d(i)*dx^3;
    k = k + 1;
end

end

