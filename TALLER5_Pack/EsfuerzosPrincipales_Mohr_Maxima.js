/**Cargamos la librería para calcular los 
 * valores y vectores propios
*/
load ("eigen");

/**Los esfuerzos son: */
sx: -1$ sy: 2$ txy: -3$
/**Armamos la matriz de tensiones */
sigma: matrix(
    [sx, txy],
    [txy, sy]
);

/**El polinomio característico se
 * calcula con:
*/
polinomcar: expand(
    charpoly(sigma, sn)
);

/**Las raíces del polinomio característico
 * son la magnitud de los esfuerzos principales
*/
Solve(
    [polinomcar=0], [sn]
);

/**Los valores y vectores propios se
 * calculan con:
*/
uniteigenvectors(sigma);
/**Y obtenemos el valor numérico */
%, numer;

/**El ángulo asociado al esfuerzo
 * principal 1 es:
*/
ang: atan2(2*txy, sx-sy)/2, numer;

/**El vector unitario asociado al esfuerzo
 * principal 1 es:
*/
[
    cos(ang), sin(ang)
], numer;
/**El vector unitario asociado al esfuerzo
 * principal 2 es:
*/
[
    cos(ang + %pi/2), sin(ang + %pi/2)
], numer;

/**Finalmente obtenemos el esfuerzo
 * cortante máximo:
*/
taumax: sqrt(
    ((sx-sy)/2)^2 + txy^2
);
%, numer;

/**El cual actúa sobre los planos cuyos
 * vectores unitarios son:
*/
[
    cos(ang + %pi/4), sin(ang + %pi/4)
], numer;
[
    cos(ang - %pi/4), sin(ang - %pi/4)
], numer;
