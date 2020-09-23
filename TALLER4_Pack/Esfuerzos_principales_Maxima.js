/**=====================================I================================== */
load ("eigen"); /**Cargamos la librería "eigen" */
sigma: matrix (
    [sx, txy, txz],
    [txy, sy, tyz],
    [txz, tyz, sz]
);

/**Definimos el vector de cosenos directores 
 * y la ecuación matricial de Cauchy
*/
ng: columnvector(
    [alpha, beta, gamma]
);
q: sigma.ng;

/**Obtenemos una expresión para el esfuerzo 
 * normal principal
*/
sigman: expand(
    transpose(q).ng
);

/**Obtenemos una expresión para el cuadrado del  
 * esfuerzo cortante principal
*/
taun2: transpose(q).q - sigman^2;

/**=====================================II================================= */
sigma: matrix (
    [sx, txy, txz],
    [txy, sy, tyz],
    [txz, tyz, sz]
);

/**Se calcula el polinomio característico y se  
 * expresa como función de sn
*/
polinomcar: expand(
    charpoly(sigma, sn)
);

/**Se extraen los coeficientes del polinomio característico, 
 * es decir, se extraen los invariantes de esfuerzo
 */
a: coeff(polinomcar, sn, 3); /**igual a (-1) */
b: coeff(polinomcar, sn, 2);
c: coeff(polinomcar, sn, 1);
d: coeff(polinomcar, sn, 0);

/**Reescribimos las líneas de código anteriores de una  
 * manera más acorde con la documentación
*/
I1: coeff(polinomcar, sn, 2);
I2: -coeff(polinomcar, sn, 1);
I3: coeff(polinomcar, sn, 0);

/**Se verifican las igualdades (2.44), por lo que cada línea 
 * debe imprimir cero (0). Aquí, el comando mat_trace(M) 
 * calcula la traza de la matriz M
 */
factor(
    I1 - mat_trace(sigma)
);
factor(
    I2 - (mat_trace(sigma)^2 - mat_trace(sigma.sigma))/2
);
factor(
    I3 - determinant(sigma)
);
