function [models, best_idx] = lsq_poly_compare(x, y, degrees, method, scaleX, doPlot)
% LSQ_POLY_COMPARE  Wrapper para comparar ajustes LS polinomicos.
%   [models,best_idx] = lsq_poly_compare(x,y,degrees,method,scaleX,doPlot)
%   x,y: datos (vectores)
%   degrees: vector de grados no negativos, ej [1 2 3 5]
%   method: 'qr' (default) o 'normal'
%   scaleX: 1 para escalar x a [-1,1], 0 para no (default: 0)
%   doPlot: 1 plotea, 0 no (default: 1)
%
%   models(i) contiene:
%     .degree  .a  .stats  .yhat  .xm  .xr
%     .xgrid   .ygrid
%
%   Ejemplo:
%     x = linspace(-2,2,25).';
%     y = sin(x) + 0.1*randn(size(x));
%     [M,idx] = lsq_poly_compare(x,y,[1 3 5],'qr',1,1);
%     disp(['Mejor grado por RMS: ', num2str(M(idx).degree)]);

if nargin < 4
    method = 'qr';
end
if nargin < 5
    scaleX = 0;
end
if nargin < 6
    doPlot = 1;
end

x = x(:);
y = y(:);

if length(x) ~= length(y)
    error('x e y deben tener la misma longitud');
end
if isempty(degrees)
    error('degrees no puede ser vacio');
end

xmin = min(x);
xmax = max(x);
if xmax == xmin
    error('Todos los x son iguales; no se puede ajustar grado > 0');
end

% Preparar grilla para curvas
xgrid = linspace(xmin, xmax, 400).';
models = struct([]);

% Parametros de escalamiento lineal a [-1,1]
xm = (xmax + xmin) / 2;
xr = (xmax - xmin) / 2;
if xr == 0
    xr = 1;
end

% Encabezado resumen
fprintf('\nResumen LS polinomico (%s). scaleX=%d\n', method, scaleX);
fprintf('Grado   SSE          RMS\n');

for k = 1:length(degrees)
    n = degrees(k);
    if n < 0
        error('Los grados deben ser >= 0');
    end

    if scaleX == 1
        xt = (x - xm) / xr;
        xgt = (xgrid - xm) / xr;
        [a, stats] = lsq_poly(xt, y, n, method);
        yhat = poly_eval_asc(a, xt);
        ygrid = poly_eval_asc(a, xgt);
        use_x = 'scaled';
    else
        [a, stats] = lsq_poly(x, y, n, method);
        yhat = poly_eval_asc(a, x);
        ygrid = poly_eval_asc(a, xgrid);
        use_x = 'raw';
    end

    models(k).degree = n;
    models(k).a = a;
    models(k).stats = stats;
    models(k).yhat = yhat;
    models(k).xgrid = xgrid;
    models(k).ygrid = ygrid;
    models(k).xm = xm;
    models(k).xr = xr;
    models(k).x_mode = use_x;

    fprintf('%5d   %11.4e   %11.4e\n', n, stats.SSE, stats.RMS);
end

% Elegir mejor por RMS
rms_vals = zeros(length(models),1);
i = 1;
while i <= length(models)
    rms_vals(i) = models(i).stats.RMS;
    i = i + 1;
end
[~, best_idx] = min(rms_vals);

% Plot
if doPlot == 1
    figure;
    plot(x, y, 'o'); hold on;
    leg = cell(length(models) + 1, 1);
    leg{1} = 'datos';
    i = 1;
    while i <= length(models)
        plot(models(i).xgrid, models(i).ygrid, '-');
        s = ['n=', num2str(models(i).degree), ' (RMS=', num2str(models(i).stats.RMS,'%.3g'), ')'];
        leg{i+1} = s;
        i = i + 1;
    end
    grid on;
    xlabel('x');
    ylabel('y');
    title('Ajustes por minimos cuadrados polinomicos');
    legend(leg, 'Location', 'Best');
    hold off;
end

end

