function labels_used = plot_xy_sets(Xcells, Ycells, connectLines, legendLocation, figTitle)
% PLOT_XY_SETS  Grafica multiples conjuntos (xi vs yi) en un solo grafico.
%   labels_used = PLOT_XY_SETS(Xcells, Ycells) grafica cada par Xcells{k}, Ycells{k}
%   como puntos (sin lineas). Xcells y Ycells deben ser cell arrays de igual
%   longitud; cada celda contiene un vector de igual tamano para ese par.
%
%   PLOT_XY_SETS(Xcells, Ycells, connectLines) si connectLines=1 une con lineas.
%   PLOT_XY_SETS(Xcells, Ycells, connectLines, legendLocation) posiciona la leyenda.
%   PLOT_XY_SETS(..., figTitle) permite definir el titulo del grafico.
%
%   Salida:
%     labels_used : cell array con los nombres efectivamente graficados.
%
%   Notas:
%   - Limpia NaN/Inf de cada par antes de graficar.
%   - Si un par tiene longitudes distintas tras limpiar, se omite y se informa.
%   - En Octave, 'best' puede no estar disponible; usar 'northeast', etc.
%
%   Ejemplo:
%     X = {x, x1, x2}; Y = {y, y1, y2};
%     plot_xy_sets(X, Y, 0, 'northeast', 'xi vs yi');
%
%   Compatibilidad: MATLAB y Octave. ASCII-only.

if nargin < 3
    connectLines = 0;
end
if nargin < 4 || isempty(legendLocation)
    legendLocation = 'northeast';
end
if nargin < 5 || isempty(figTitle)
    figTitle = 'Conjuntos xi vs yi';
end

if ~iscell(Xcells) || ~iscell(Ycells)
    error('Xcells y Ycells deben ser cell arrays.');
end
if numel(Xcells) ~= numel(Ycells)
    error('Xcells y Ycells deben tener igual cantidad de celdas.');
end

mk = {'o','s','^','d','v','>','<','x','+','p','h'};

figure;
hold on;
labels_used = {};
idxlab = 1;

k = 1;
while k <= numel(Xcells)
    xi = Xcells{k};
    yi = Ycells{k};

    if isempty(xi) || isempty(yi)
        fprintf('Conjunto %d vacio -> skip\n', k-1);
        k = k + 1;
        continue;
    end

    xi = xi(:);
    yi = yi(:);

    mask = isfinite(xi) & isfinite(yi);
    xi = xi(mask);
    yi = yi(mask);

    if numel(xi) ~= numel(yi)
        fprintf('Conjunto %d mismatch tras limpiar: %d vs %d -> skip\n', k-1, numel(xi), numel(yi));
        k = k + 1;
        continue;
    end
    if numel(xi) == 0
        fprintf('Conjunto %d sin datos finitos -> skip\n', k-1);
        k = k + 1;
        continue;
    end

    mki = mk{mod(k-1, numel(mk)) + 1};
    if connectLines == 1
        plot(xi, yi, mki, 'MarkerSize', 6, 'LineStyle', '-');
    else
        plot(xi, yi, mki, 'MarkerSize', 6, 'LineStyle', 'none');
    end

    labels_used{idxlab,1} = ['Conjunto ', num2str(k-1)];
    idxlab = idxlab + 1;

    k = k + 1;
end

grid on;
xlabel('x');
ylabel('y');
title(figTitle);

if ~isempty(labels_used)
    legend(labels_used, 'Location', legendLocation);
end

hold off;

end
