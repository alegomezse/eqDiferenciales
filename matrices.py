import numpy.matlib 
import numpy as np 
from scipy.linalg import lu
from fractions import Fraction

def raicesPolinomio(a,b,c,d):
    divisores = []
    raices = []
    for i in range(abs(d)+1):
        if(i == 0):
            pass
        else:
            if(d % i==0):
                divisores.append(i)
    for j in range(len(divisores)):
        divisores.append(-divisores[j])

    print("divisores = ",divisores)
    for k in range(len(divisores)):
        var1 = var2 = var3 = var4 = 0
        var1 = b + (divisores[k]*a)
        var2 = divisores[k]*var1
        var3 = var2 + c
        var4 = divisores[k]*var3 + d
        if(var4 == 0):
            raices.append(divisores[k])
    print("raices = ",raices)

def ecuacionValoresPropiosR3(matriz):
    #dada una matriz 3x3 se calula la ecuacion cubica correspondiente para los valores propios
    t = matriz[0,0] + matriz[1,1] + matriz[2,2]
    B = (matriz[0,0]*matriz[1,1]+matriz[0,0]*matriz[2,2]+matriz[1,1]*matriz[2,2])-(matriz[0,1]*matriz[1,0]+matriz[0,2]*matriz[2,0]+matriz[1,2]*matriz[2,1])
    g = matriz[0,0]*(matriz[1,1]*matriz[2,2]-matriz[1,2]*matriz[2,1])-matriz[0,1]*(matriz[1,0]*matriz[2,2]-matriz[1,2]*matriz[2,0])+matriz[0,2]*(matriz[1,0]*matriz[2,1]-matriz[2,0]*matriz[1,1])
    ecuacionCubica = "ecuacion = y^3 - "+str(Fraction(t).limit_denominator())+"y^2  + ("+str(Fraction(B).limit_denominator())+")y -("+str(Fraction(g).limit_denominator())+")"
    print(ecuacionCubica)
    raicesPolinomio(1,-int(round(t)),int(round(B)),-int(round(g)))
    
def calculoVectorImpropioR3(matriz,v1):
    #dado un valor impropio y la  matriz original se calculan la matriz reducida para encontrar los valores
    matrizv1 = np.array([[v1,0,0],[0,v1,0],[0,0,v1]])
    matrizValorImpropio = matriz - matrizv1
    print("matriz valor impropio valor "+str(v1)+" = ")
    print(matrizValorImpropio)
    pl, matrizGaussJordan = lu(matrizValorImpropio, permute_l=True)
    print("matriz gauss-jordan valor "+str(v1)+" = ")
    print(matrizGaussJordan)
    print("")
def vectorAsociadoR3(matriz,v1):
    identidad = np.array([[1,0,0],[0,1,0],[0,0,1]])
    ecuacion1 = matriz - identidad
    ecuacion2 = np.insert(ecuacion1, ecuacion1.shape[1],np.array(v1), 1)
    pl, MgaussJordan = lu(ecuacion2, permute_l=True)
    print("gauss jordan vector asociado ")
    print(MgaussJordan)
def matrizJCJ(matriz,matrizVectoresImpropiosR3):
    inverse = np.linalg.inv(matrizVectoresImpropiosR3)
    jc = np.dot(inverse,matriz)
    jcj = np.dot(jc,matrizVectoresImpropiosR3)  
    print("matriz j = ")
    print(jcj)


def matrizR2(matriz2x2):
   #dada una matriz 2x2 se calculan sus valor propios y la matriz reducida respectivamente
   B = -1*(matriz2x2[0,0]+matriz2x2[1,1])
   C = matriz2x2[0,0]*matriz2x2[1,1]-matriz2x2[0,1]*matriz2x2[1,0]
   vPositivo = (-B + np.sqrt(pow(B,2)-4*C))/2
   vNegativo = (-B - np.sqrt(pow(B,2)-4*C))/2
   print("valor positivo = ",vPositivo)
   print("valor negativo = ",vNegativo)
   matrizValorImpropioPositivo = np.array([ [matriz2x2[0,0]-vPositivo,matriz2x2[0,1]] , [matriz2x2[1,0],matriz2x2[1,1]-vPositivo]])
   matrizValorImpropioNegativo = np.array([ [matriz2x2[0,0]-vNegativo,matriz2x2[0,1]] , [matriz2x2[1,0],matriz2x2[1,1]-vNegativo]])
   #print(matrizValorImpropioPositivo)
   #print(matrizValorImpropioNegativo)
   pl, gaussValorpositivo = lu(matrizValorImpropioPositivo, permute_l=True)
   pl, gaussValorNegativo = lu(matrizValorImpropioNegativo, permute_l=True) 
   print("ecuacion gauss jordan 1 :")
   print(gaussValorpositivo)
   print("ecuacion gauss jordan 2 :")
   print(gaussValorNegativo)
 
def solucionGeneralR2(matrizVectoresImpropios):
    #dada la matriz de vectores propios se calcula la ecuación con la solución general
    inverse = np.linalg.inv(matrizVectoresImpropios)
    matrizInversaTexto = [[str(Fraction(inverse[0,0]).limit_denominator())+"k1 "+"+ ("+str(Fraction(inverse[0,1]).limit_denominator())+"k2)"],[str(Fraction(inverse[1,0]).limit_denominator())+"k1 "+" + ("+str(Fraction(inverse[1,1]).limit_denominator())+"k2)"]]
    matrizTexto = [[str(matrizVectoresImpropios[0,0]),str(matrizVectoresImpropios[0,1])],[str(matrizVectoresImpropios[1,0]),str(matrizVectoresImpropios[1,1])]]
    arrayEuler = ([[matrizTexto[0][0]+"e^(v1*t)",matrizTexto[0][1]+"e^v2*t"],[matrizTexto[1][0]+"e^(v1*t)",matrizTexto[1][1]+"e^v2*t"]])
    print("")
    print("a :",arrayEuler)
    print("")
    print("b :",matrizInversaTexto)
    print("")
    print("multiplique a * b")

matriz2x2 = np.array([[3.0,-2.0],[2.0,-2.0]])#matriz A en R2
matrizVectoresImpropios = np.array([[2.0,1.0],[1.0,2.0]])#Matriz de vectores impropios

matriz = np.array([[1.0,0.0,0.0],[-4.0,1.0,0.0],[3.0,6.0,2.0]])#matriz A en R3
matrizVectoresImpropiosR3 = np.array([[0,-1/4,0],[1,-7/8,0],[-6,0,1]])#matriz vectores impropios
vectorA = np.array([0,1,-6])
ecuacionValoresPropiosR3(matriz)#calcula los valores propios y la ecuacion
calculoVectorImpropioR3(matriz,1)#encuentra la matriz gauss jordan del vector propio respectivo   
calculoVectorImpropioR3(matriz,2) 

vectorAsociadoR3(matriz,vectorA)
matrizJCJ(matriz,matrizVectoresImpropiosR3)

#matrizVectoresImpropiosr3 = np.array([[1,1,2],[-4,0,1],[1,-1,2]])
#raicesPolinomio(1,-6,-15,-8)

 
#matrizR2(matriz2x2)
#solucionGeneralR2(matrizVectoresImpropios)
#ecuacionValoresPropiosR3(matriz)
#calculoVectorImpropioR3(matriz,v1)
#inversaR3(matrizVectoresImpropiosR3)