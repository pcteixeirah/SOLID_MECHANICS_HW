import numpy as np
'''
Considere un punto P sometido a los esfuerzos: sigma_x, sigma_y,
sigma_z, tao_xy, tao_xz, tao_yz
   ---> Encuentre las direcciones, las magnitudes y los planos
   sobre los que actúan los esfuerzos principales en dicho punto
'''
# =============================================================================
def null(A, eps=1e-15):
    u, s, vh = np.linalg.svd(A)
    null_space = np.compress(s <= eps, vh, axis=0)
    return null_space.T
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
result2 = null(A2).T
# =============================================================================
# #Esfuerzo principal sigma_3
# =============================================================================
sigma_3 = poly_roots[2]
A3 = sigma_matrix
A3 = A3-sigma_3*np.identity(3)
result3 = solution(A3)
# =============================================================================
# RESULTADOS
# =============================================================================
resulti = np.cross(result1, result2)
print('¿Conforman los vectores propios una base ortogonal?: ')
print(result3)
print(resulti)
print('¿Son iguales estos dos vectores?: ')
query=eval(input("---------(¿1/0?)--------->"))
if query>0:
    eugenval = np.identity(3)*poly_roots
    print("Esfuerzos principales: \n", eugenval," [Pa]")
    print("Direcciones principales: \n", result1, " (Principal mayor) \n",
          result2, "(Principal menor) \n", result3, "(Negativo) \n")
else:
    print('FATAL ERROR')









