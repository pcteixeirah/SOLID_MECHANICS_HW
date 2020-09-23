import numpy as np
'''
Considere un punto de un sólido bidimensional en el cual los esfuerzos
en un punto dado son: sigma_x, sigma_y, tao_xy
   1--Plantear la matriz de tensiones sigma_matrix correspondiente
   2--Calcular el polinomio característico asociado a sigma_matrix
   3--Calcular la dirección y magnitud de los esfuerzos principales
'''
sigma_x = 3
sigma_y = 2
tao_xy = -4

#Definimos la matriz de tensiones sigma_matrix
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
eugenval = np.zeros((2,2))
eugenval[0,0] = sigma_min
eugenval[1,1] = sigma_max

eugenvect = np.transpose(np.mat([n2_eugenvect, n1_eugenvect]))

print("Esfuerzos principales: \n", eugenval," [Pa]")
print("Direcciones principales: \n", eugenvect)

