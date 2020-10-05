import numpy as np
'''
Considere un punto P sometido a los esfuerzos: sigma_x, sigma_y,
sigma_z, tao_xy, tao_xz, tao_yz
   ---> Encuentre las direcciones y los esfuerzos principales en dicho punto
   ---> Encuentre los vectores normalizados que definen los planos para los
   cuales se presenta el estado de esfuerzos: sigma_n, tao_n
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
#Determinamos el estado de esfuerzos requerido en el plano v:
sigma_n = 0.1
tao_n = 2
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
#Calculamos el valor absoluto de los cosenos directores correspondientes a los
# planos v según las ecuaciones 2.68
arg_1 = (tao_n**2+(sigma_n-sigma_2)*(sigma_n-sigma_3))/((sigma_1-sigma_2)*(sigma_1-sigma_3))
alpha = np.abs(np.sqrt(arg_1))

arg_2 = (tao_n**2+(sigma_n-sigma_3)*(sigma_n-sigma_1))/((sigma_2-sigma_3)*(sigma_2-sigma_1))
beta = np.abs(np.sqrt(arg_2))

arg_3 = (tao_n**2+(sigma_n-sigma_1)*(sigma_n-sigma_2))/((sigma_3-sigma_1)*(sigma_3-sigma_2))
gamma = np.abs(np.sqrt(arg_3))

#Definimos los 4 vectores asociados a las combinaciones de los cosenos 
# directores calculados anteriormente, recuérdese que dichos vectores
# están definidos con respecto a la base ortogonal de las direcciones principales
v1 = np.transpose([alpha, beta, gamma])
v2 = np.transpose([alpha, beta, -gamma])
v3 = np.transpose([alpha, -beta, gamma])
v4 = np.transpose([alpha, -beta, -gamma])
# =============================================================================
print("Esfuerzos principales: \n", poly_roots," [Pa]")
print("Direcciones principales: \n", result1, " (Principal mayor) \n",
          result2, "(Principal media) \n", result3, "(Principal menor) \n")
print("Dirección 1 correspondiente a sigma_n, tao_n: \n", v1)
print("Dirección 2 correspondiente a sigma_n, tao_n: \n", v2)
print("Dirección 3 correspondiente a sigma_n, tao_n: \n", v3)
print("Dirección 4 correspondiente a sigma_n, tao_n: \n", v4)

