import numpy as np
import matplotlib.pyplot as plt
'''
Considere un punto P sometido a los esfuerzos: sigma_x, sigma_y,
sigma_z, tao_xy, tao_xz, tao_yz
   ---> Encuentre las direcciones y los esfuerzos principales en dicho punto
   ---> Encuentre los vectores normalizados que definen los planos en que se
   presentan los cortantes máximos, así como la magnitud de dichos esfuerzos
   ---> Grafique el diagrama de círculos de Mohr correspondiente a este 
   estado de esfuerzos
'''
# =============================================================================
def solution(U):
    #Find the eigenvalues and eigenvector of U(transpose).U
    e_vals, e_vecs = np.linalg.eig(np.dot(U.T, U))  
    #Extract the eigenvector (column) associated with the 
    # minimum eigenvalue
    return e_vecs[:, np.argmin(e_vals)] 
# =============================================================================
sigma_x = 1
sigma_y = 3
sigma_z = 0
tao_xy = 0
tao_xz = 0
tao_yz = 2
# =============================================================================
#Definimos la matriz de tensiones sigma_matrix
sigma_matrix = np.zeros((3,3))
sigma_matrix[0,0] = sigma_x
sigma_matrix[0,1] = tao_xy
sigma_matrix[0,2] = tao_xz
sigma_matrix[1,0] = tao_xy
sigma_matrix[1,1] = sigma_y
sigma_matrix[1,2] = tao_yz
sigma_matrix[2,0] = tao_xz
sigma_matrix[2,1] = tao_yz
sigma_matrix[2,2] = sigma_z

#Definimos los coeficientes del polinomio característico
# asociado a sigma_matrix
# a*sigma_n**3+b*sigma_n**2+c*sigma_n+d = 0
char_poly = np.poly(sigma_matrix)

#Definimos las raíces del polinomio en sigma_n
poly_roots = np.roots(char_poly)
poly_roots = np.sort(poly_roots)[::-1] #Ordenamos de mayor a menor

#Resolvemos los sistemas de ecuaciones a fin de obtener los
# vectores asociados a las direcciones principales
# =============================================================================
# #Esfuerzo principal sigma_1
# =============================================================================
sigma_1 = poly_roots[0]
A1 = sigma_matrix
A1 = A1-sigma_1*np.identity(3)
result1 = solution(A1)
# =============================================================================
# #Esfuerzo principal sigma_2
# =============================================================================
sigma_2 = poly_roots[1]
A2 = sigma_matrix
A2 = A2-sigma_2*np.identity(3)
result2 = solution(A2)
# =============================================================================
# #Esfuerzo principal sigma_3
# =============================================================================
sigma_3 = poly_roots[2]
A3 = sigma_matrix
A3 = A3-sigma_3*np.identity(3)
result3 = solution(A3)
# =============================================================================
#Determinamos el esfuerzo cortante máximo y los vectores que definen los 
# planos en los que se presenta, según las expresiones conocidas 
tao_max = (sigma_1-sigma_3)/2
v1_diag = (result1+result3)/np.linalg.norm(result1+result3)
v2_diag = (result1-result3)/np.linalg.norm(result1-result3)

# =============================================================================
# DIAGRAMA DE CÍRCULOS DE MOHR 3D
# =============================================================================

# =============================================================================
#Definimos los parámetros de la circunferencia exterior C2
delta = 0.005 #Accuracy
r2 = (sigma_1-sigma_3)/2   #Radio de la semicircunferencia MOHR
h2 = (sigma_1+sigma_3)/2 #Centro de la semicircunferencia MOHR

sigma_mohr2 = np.arange(sigma_3, sigma_1+delta, delta) #Esfuerzos normales MOHR
n2_parameter = np.size(sigma_mohr2)
ang2_parameter = np.linspace(0, 2*np.pi, n2_parameter+1)

sigma_c2 = r2*np.cos(ang2_parameter)+h2
tao_c2 = r2*np.sin(ang2_parameter)
# =============================================================================
#Definimos los parámetros de la circunferencia interior C1
r1 = (sigma_2-sigma_3)/2   #Radio de la semicircunferencia MOHR
h1 = (sigma_2+sigma_3)/2 #Centro de la semicircunferencia MOHR

sigma_mohr1 = np.arange(sigma_3, sigma_2+delta, delta) #Esfuerzos normales MOHR
n1_parameter = np.size(sigma_mohr1)
ang1_parameter = np.linspace(0, 2*np.pi, n1_parameter+1)

sigma_c1 = r1*np.cos(ang1_parameter)+h1
tao_c1 = r1*np.sin(ang1_parameter)
# =============================================================================
#Definimos los parámetros de la circunferencia interior C3
r3 = (sigma_1-sigma_2)/2   #Radio de la semicircunferencia MOHR
h3 = (sigma_1+sigma_2)/2 #Centro de la semicircunferencia MOHR

sigma_mohr3 = np.arange(sigma_2, sigma_1+delta, delta) #Esfuerzos normales MOHR
n3_parameter = np.size(sigma_mohr3)
ang3_parameter = np.linspace(0, 2*np.pi, n3_parameter+1)

sigma_c3 = r3*np.cos(ang3_parameter)+h3
tao_c3 = r3*np.sin(ang3_parameter)
# =============================================================================
#Definimos los puntos en la gráfica que corresponden a los esfuerzos 
# principales (normales) y cortantes máximos
point_sigma1 = [sigma_1, 0]
point_sigma2 = [sigma_2, 0]
point_sigma3 = [sigma_3, 0]
point_taomax = [h2, tao_max]
point_taomin = [h2, -tao_max]
# =============================================================================
print("Esfuerzos principales: \n", poly_roots," [Pa]")
print("Esfuerzo cortante máximo: \n", tao_max," [Pa]")
print("Direcciones principales: \n", result1, " (Principal mayor) \n",
          result2, "(Principal media) \n", result3, "(Principal menor) \n")
print("Dirección 1 correspondiente a tao_max: \n", v1_diag)
print("Dirección 2 correspondiente a tao_max: \n", v2_diag)
# =============================================================================
plt.axis('scaled')
plt.ylim(-r2-0.3, r2+0.7) 
plt.xlim(sigma_3-0.3, sigma_1+0.7)
plt.grid()

plt.plot(sigma_c2, tao_c2)
plt.plot(sigma_c1, tao_c1)
plt.plot(sigma_c3, tao_c3)

plt.plot(point_sigma1[0],point_sigma1[1],'go')
plt.annotate('S1', point_sigma1, xytext=(4,6),textcoords='offset pixels')
plt.plot(point_sigma2[0],point_sigma2[1],'go')
plt.annotate('S2', point_sigma2, xytext=(4,6),textcoords='offset pixels')
plt.plot(point_sigma3[0],point_sigma3[1],'go')
plt.annotate('S3', point_sigma3, xytext=(4,6),textcoords='offset pixels')
plt.plot(point_taomax[0],point_taomax[1],'go')
plt.annotate('Tao', point_taomax, xytext=(4,6),textcoords='offset pixels')
plt.plot(point_taomin[0],point_taomin[1],'go')
plt.annotate('Tao', point_taomin, xytext=(4,6),textcoords='offset pixels')

plt.show()
# =============================================================================


