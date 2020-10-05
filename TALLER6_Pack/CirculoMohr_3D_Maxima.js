/**========================================I========================================*/
/**Definimos nuestro sistema de ecuaciones,
 * compuesto por 2.65, 2.66, 2.67
*/
eq1: sn = s1*a^2 + s2*b^2 + s3*g^2$
eq2: tn^2 = (s1*a)^2 + (s2*b)^2 + (s3*g)^2 - sn^2$
eq3: a^2 + b^2 + g^2 = 1$

/**Resolvemos para alpha, beta y gamma */
solve(
    [eq1, eq2, eq3],
    [a, b, g]
)$
factor(%^2)$ /**Se calculan y factorizan los cuadrados*/

/**Aquí aparecerán varios resultados en listas
 * diferentes, algunos de ellos repetidos, siendo
 * necesario quitarlos
*/
flatten(%)$ /**Se ponen los resultados en una lista */
transpose(
    unique(%)
); /**Se quitan las soluciones repetidas */

/**========================================II=======================================*/
/**Se definen los esfuerzos normales y cortantes en el
 * punto de cortante máximo, en función de los esfuerzos
 * normales principales (s1, s2, s3)
*/
sn: (s1 + s3)/2$
tn: (s1 - s3)/2$

/**Conociendo las expresiones para los cosenos directores
 * alpha, beta y gamma, dadas por 2.68, podemos calcular
 * las componentes del vector normal asociado a las 
 * condiciones de cortante máximo
*/
alpha2: factor((tn^2+sn^2-s3*sn-s2*sn+s2*s3)/((s2-s1)*(s3-s1)));
beta2: factor(-(tn^2+sn^2-s3*sn-s1*sn+s1*s3)/((s2-s1)*(s3-s2)));
gamma2: factor((tn^2+sn^2-s2*sn-s1*sn+s1*s2)/((s3-s1)*(s3-s2)));