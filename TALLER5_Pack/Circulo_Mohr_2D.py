import numpy as np
import matplotlib.pyplot as plt
import principal_directions as pr
# =============================================================================
'''
Considere un punto sujeto a los esfuerzos:
sigma_x, sigma_y, tao_xy [Pa]
Encuentre los esfuerzos principales (y su dirección) 
para el punto en consideración
'''
# =============================================================================
sigma_x = -1
sigma_y = 2
tao_xy = -3
# =============================================================================
delta = 0.005 #Accuracy
r = np.sqrt(((sigma_x-sigma_y)/2)**2+tao_xy**2)   #Radio de la semicircunferencia MOHR
h = (sigma_x+sigma_y)/2 #Centro de la semicircunferencia MOHR
sigma_min = h-r
sigma_max = h+r
sigma_mohr = np.arange(sigma_min, sigma_max+delta, delta) #Esfuerzos normales MOHR
n_parameter = np.size(sigma_mohr)
ang_parameter = np.linspace(0, 2*np.pi, n_parameter+1);

sigma_axis = r*np.cos(ang_parameter)+h;
tao_axis = r*np.sin(ang_parameter);

mohr_circle = list(zip(sigma_axis, tao_axis))

point_hk = [h, 0]
point_A = [sigma_x, tao_xy]
point_C = [sigma_y, -tao_xy]
point_max = [sigma_max, 0]
point_min = [sigma_min, 0]

# =============================================================================
oplane_sigma = [point_A[0],point_C[0]]
oplane_tao = [point_A[1],point_C[1]]

op_vector = np.subtract(point_A,point_C)
norm_op = np.linalg.norm(op_vector)
unitvec_op = np.divide(op_vector,norm_op)
# =============================================================================
flatplane_sigma = [point_max[0],point_min[0]]
flatplane_tao = [point_max[1],point_min[1]]

flat_vector = np.subtract(point_max,point_min)
norm_flat = np.linalg.norm(flat_vector)
unitvec_flat = np.divide(flat_vector,norm_flat)
# =============================================================================
theta1_double = np.arctan2(tao_xy,sigma_x-h)
theta1_double = np.rad2deg(theta1_double)
theta1 = theta1_double/2
theta2 = theta1+90
thetac1 = theta1-45
thetac2 = theta1+45
# =============================================================================
#Definimos la matriz de tensiones sigma_matrix
sigma_matrix = np.zeros((2,2))
sigma_matrix[0,0] = sigma_x
sigma_matrix[0,1] = tao_xy
sigma_matrix[1,0] = tao_xy
sigma_matrix[1,1] = sigma_y

#Definimos los coeficientes del polinomio característico
char_poly = np.poly(sigma_matrix)
#Definimos las raíces del polinomio en sigma_n
poly_roots = np.roots(char_poly)
poly_roots = np.sort(poly_roots)[::-1] #Ordenamos de mayor a menor

eugenvect = pr.dir2D(sigma_x,sigma_y,tao_xy)

identity = np.identity(2)
sigma_index = poly_roots[0]
sigma_identity = sigma_index*identity
sigma_term = np.subtract(sigma_matrix, sigma_identity)
ng_vect = np.zeros((2,1))
ng_vect[0,0] = eugenvect[0,0]
ng_vect[1,0] = eugenvect[1,0]
tao_max = np.matmul(sigma_term, ng_vect)
tao_max = np.linalg.norm(tao_max)/2

theta_num = np.arccos(eugenvect[0,0])
theta_num = np.rad2deg(theta_num)
if tao_xy<0:
    theta2_num = theta_num
    theta1_num = theta2_num-90
    thetac1_num = theta1_num-45
    thetac2_num = theta1_num+45
else:
    theta1_num = theta_num
    theta2_num = theta1_num+90
    thetac1_num = theta1_num-45
    thetac2_num = theta1_num+45
# =============================================================================
print("Esfuerzos principales (numérico): \n", poly_roots," [Pa]")
print("Esfuerzos principales (gráfico): \n [", 
      point_max[0], point_min[0], "] [Pa]")
print("Esfuerzo cortante (numérico): \n", tao_max," [Pa]")
print("Esfuerzo cortante (gráfico): \n", r," [Pa]")
print("Direcciones principales: \n", eugenvect)
print("Theta_1 [deg] (numérico): \n", theta1_num)
print("Theta_1 [deg] (gráfico): \n", theta1)
print("Theta_2 [deg] (numérico): \n", theta2_num)
print("Theta_2 [deg] (gráfico): \n", theta2)
print("Theta_c1 [deg] (numérico): \n", thetac1_num)
print("Theta_c1 [deg] (gráfico): \n", thetac1)
print("Theta_c2 [deg] (numérico): \n", thetac2_num)
print("Theta_c2 [deg] (gráfico): \n", thetac2)
# =============================================================================
plt.axis('scaled')
plt.ylim(-r-0.3, r+0.7) 
plt.xlim(sigma_min-0.3, sigma_max+0.7)
plt.grid()

plt.plot(sigma_axis, tao_axis)

plt.plot(point_A[0],point_A[1],'go')
plt.annotate('A', point_A, xytext=(4,6),textcoords='offset pixels')
plt.plot(point_C[0],point_C[1],'go')
plt.annotate('C', point_C, xytext=(4,6),textcoords='offset pixels')
plt.plot(point_hk[0],point_hk[1],'go')
plt.plot(point_max[0],point_max[1],'go')
plt.annotate('B', point_max, xytext=(4,6),textcoords='offset pixels')
plt.plot(point_min[0],point_min[1],'go')
plt.annotate('D', point_min, xytext=(4,6),textcoords='offset pixels')

plt.plot(oplane_sigma,oplane_tao)
plt.plot(flatplane_sigma,flatplane_tao)

plt.show()
# =============================================================================
