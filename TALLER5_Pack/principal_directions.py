import numpy as np
# =============================================================================
def dir2D(sigma_x, sigma_y, tao_xy):
    sigma_matrix = np.zeros((2,2))
    sigma_matrix[0,0] = sigma_x
    sigma_matrix[0,1] = tao_xy
    sigma_matrix[1,0] = tao_xy
    sigma_matrix[1,1] = sigma_y

#Definimos los coeficientes del polinomio característico
# asociado a sigma_matrix
# a*sigma_n**2+b*sigma_n+c = 0
    a = 1
    b = -(sigma_x+sigma_y)
    c = sigma_x*sigma_y-tao_xy**2

#Definimos las raíces del polinomio en sigma_n
    sigma_1 = (-b+np.sqrt(b**2-4*a*c))/(2*a)
    sigma_2 = (-b-np.sqrt(b**2-4*a*c))/(2*a)
    sigma_max = np.max([sigma_1, sigma_2])
    sigma_min = np.min([sigma_1, sigma_2])

#Resolvemos el sistema de ecuaciones 2.40 a fin de obtener los
# vectores asociados a las direcciones principales
# =============================================================================
# #Esfuerzo principal mayor
# =============================================================================
    alpha1_coef = (sigma_x-sigma_max)
    beta1_coef = tao_xy
#Resolvemos para beta1 = 1
    beta1 = 1
    alpha1 = -beta1_coef/alpha1_coef
    n1 = np.transpose([alpha1, beta1])
#Normalizamos para obtener el vector propio
    n1_eugenvect = np.divide(n1, np.linalg.norm(n1))
# =============================================================================
# #Esfuerzo principal menor
# =============================================================================
    alpha2_coef = (sigma_x-sigma_min)
    beta2_coef = tao_xy
#Resolvemos para beta2 = 1
    beta2 = 1
    alpha2 = -beta2_coef/alpha2_coef
    n2 = np.transpose([alpha2, beta2])
#Normalizamos para obtener el vector propio
    n2_eugenvect = np.divide(n2, np.linalg.norm(n2))

# =============================================================================
# RESULTADOS
# =============================================================================
    eugenvect = np.transpose(np.mat([n2_eugenvect, n1_eugenvect]))
    return eugenvect

# =============================================================================
def solution(U):
    #Find the eigenvalues and eigenvector of U(transpose).U
    e_vals, e_vecs = np.linalg.eig(np.dot(U.T, U))  
    #Extract the eigenvector (column) associated with the 
    # minimum eigenvalue
    return e_vecs[:, np.argmin(e_vals)] 
# =============================================================================


