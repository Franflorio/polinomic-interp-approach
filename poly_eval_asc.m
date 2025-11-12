function y = poly_eval_asc(a, x)
% POLY_EVAL_ASC  Evalua p(x)=a0 + a1 x + ... + an x^n con Horner ascendente.
x = x(:);
n = length(a) - 1;
y = a(n+1) * ones(size(x));
k = n;
while k >= 1
    y = a(k) + x .* y;
    k = k - 1;
end
end

