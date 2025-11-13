# Unidad 5 — Interpolación y Aproximación de Funciones (MATLAB/Octave)

Este repositorio reúne implementaciones claras (compatibles MATLAB/Octave, solo ASCII) de:

- **Interpolación polinómica:** Lagrange (clásico y baricéntrico), Newton (diferencias divididas), Coeficientes indeterminados (Vandermonde).
- **Splines cúbicos:** natural y sujeto (clamped).
- **Aproximación por mínimos cuadrados:** ecuaciones normales y variante estable con `\` (QR).

> **Fuentes**: teoría y fórmulas alineadas con los **apuntes de cátedra** que subiste. La función `chebyshev_nodes.m` se incluye como **utilidad opcional** estándar de la bibliografía, **no aparece** explícitamente en los PDFs de la cátedra; podés omitirla si tu programa no la requiere.

---

## Requisitos y convenciones

- MATLAB **o** Octave.
- Solo ASCII, sin operadores ternarios ni continuaciones de línea `...`.
- Vectores de entrada aceptan fila o columna; internamente se normalizan a **columna**.
- Polinomios en base monómica usan **coeficientes ascendentes**: `p(x) = a0 + a1 x + ... + an x^n`.

---

## Estructura sugerida

```
/src
  barycentric_eval.m
  barycentric_weights.m
  chebyshev_nodes.m
  interp_eval_point.m
  lagrange_eval_basic.m
  newton_dd.m
  newton_eval.m
  poly_eval_asc.m
  spline_cubica.m
  spline_eval.m
  vandermonde_interp.m
  lsq_poly_normales.m
  lsq_poly.m
  lsq_poly_compare.m
  plot_xy_sets.m

/examples
  u5_demo.m
  tp05_lsq_driver.m
```

> Si no usás todas las utilidades, dejá solo las necesarias para tu TP.

---

## Métodos y archivos (resumen)

### Lagrange
- **Clásico:** `lagrange_eval_basic.m` — didáctico, costo O(m n^2).
- **Baricéntrico:** `barycentric_weights.m` + `barycentric_eval.m` — estable y O(m n).

**Uso**
```matlab
w = barycentric_weights(x);
p = barycentric_eval(x, y, w, xq);
```

### Newton (diferencias divididas)
- `newton_dd.m` (coeficientes) y `newton_eval.m` (evaluación anidada).  
**Ventaja:** fácil de actualizar si agregás un punto.

### Coeficientes indeterminados (Vandermonde)
- `vandermonde_interp.m` — resuelve `V*a = y`. Único interpolante si `x` son distintos (puede ser mal condicionada para grados altos).

### Splines cúbicos
- `spline_cubica.m` (natural o clamped) + `spline_eval.m`.  
**Idea:** cúbicos por tramo, continuidad de S, S', S''; sistema tridiagonal en las curvaturas.

### Mínimos cuadrados polinomiales
- `lsq_poly_normales.m` — ecuaciones normales (reporta SSE, RMS y sigma^2).
- `lsq_poly.m` — variante estable con `\` (QR/backslash).
- `lsq_poly_compare.m` — recorre grados, opcionalmente escala x, plotea y elige por RMS.
- `tp05_lsq_driver.m` — recorre grados 0:nmax, elige **menor varianza** (sigma^2) y grafica.

### Utilidades
- `poly_eval_asc.m` — Horner ascendente (rápido).
- `chebyshev_nodes.m` — **opcional** (no figura en tus PDFs); útil para mitigar Runge en interpolación global.
- `interp_eval_point.m` — evalúa en **un punto** (base Newton o monomial).
- `plot_xy_sets.m` — ploteo robusto de muchos conjuntos (xi vs yi).

---

## Recetas rápidas

**Interpolar con baricéntrico**
```matlab
w = barycentric_weights(x);
xq = linspace(min(x), max(x), 400).';
p  = barycentric_eval(x, y, w, xq);
plot(xq, p, '-', x, y, 'o', 'LineStyle', 'none'), grid on
```

**Interpolar con Newton**
```matlab
[aN,~] = newton_dd(x, y);
pN = newton_eval(x, aN, xq);
```

**Spline natural**
```matlab
S  = spline_cubica(x, y, 'natural');
ps = spline_eval(S, xq);
```

**Mínimos cuadrados (grado 3)**
```matlab
[a3,s3] = lsq_poly_normales(x, y, 3);
p3 = poly_eval_asc(a3, xq);
```

**Driver del TP**
```matlab
run tp05_lsq_driver.m
```

---

## Checklist

- [x] Lagrange (clásico y baricéntrico)
- [x] Newton (diferencias divididas)
- [x] Vandermonde
- [x] Splines (natural / clamped)
- [x] Mínimos cuadrados (normales y QR)
- [x] Utilidades (Horner, plotting, etc.)
- [ ] Opcional: LS ponderado / ridge
- [ ] Opcional: derivadas del spline `S'`, `S''`
- [ ] Opcional: Newton incremental (agregar puntos sin recomputar)

---

## Licencia y créditos

- Código educativo. Teoría alineada con tus apuntes de cátedra.
- `chebyshev_nodes.m` es una utilidad estándar de la bibliografía, no obligatoria según tu programa.
