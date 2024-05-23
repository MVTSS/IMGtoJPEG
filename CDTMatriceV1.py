"""
Réduction d'une matrice avec la transformée en cosinus discret
Matisse A.

v1
"""

import math as ma

#Matrice Base (donnée dans le TD)
Mbase = [[70, 70, 100, 70, 87, 87, 150, 187],
         [85, 100, 96, 79, 87, 154, 87, 113],
         [100, 85, 116, 79, 70, 87, 86, 196],
         [136, 69, 87, 200, 79, 71, 117, 96],
         [161, 70, 87, 200, 103, 71, 96, 113],
         [161, 123, 147, 133, 113, 113, 85, 161],
         [146, 147, 175, 100, 103, 103, 163, 187],
         [156, 146, 186, 70, 113, 161, 163, 197]]

# Résultat voulu pour comparaison après CDT
ResVoulu = [[-80, -40, 89, -73, 44, 33, 53, -3],
            [-135, -59, -26, 6, 14, -3, -13, -28],
            [47, -76, 66, -3, -108, -78, 33, 59],
            [-2, 10, -18, 0, 33, 11, -21, 1],
            [-1, -9, -22, 8, 32, 65, -36, -1],
            [5, -20, 28, -46, 3, 24, -30, 24],
            [6, -20, 37, -28, 12, -35, 33, 17],
            [-5, -23, 33, -30, 17, -5, -4, 20]]

# Tableau du PasQuantitatif pour Cq
PasQuantif = [[16, 11, 10, 16, 24, 40, 51, 61],
              [12, 12, 14, 19, 26, 58, 60, 55],
              [14, 13, 16, 24, 40, 57, 69, 56],
              [14, 17, 22, 29, 51, 87, 80, 62],
              [18, 22, 37, 56, 69, 109, 103, 77],
              [24, 35, 55, 64, 81, 104, 113, 92],
              [49, 64, 78, 87, 103, 121, 120, 101],
              [72, 92, 95, 98, 112, 100, 103, 99]]




def f(u,v, base):
    """Prend l'élément (u,v) d'une matrice "base"

    Args:
        u (int): équivalent de x
        v (int): équivalent de y
        base (list): matrice carré

    Returns:
        int: élément correspondant à l'emplacement (u,v) de la matrice
    """
    return base[u][v]


def c(k):
    if k == 0 :
        return 1/ma.sqrt(2) 
    else:
        return 1

def CDT(u,v, matrice):
    fcos = 0
    x, y, N = c(u), c(v), len(matrice)
    for n in range(N):
        for m in range(N):
            fcos += f(n,m, matrice)*ma.cos(((2*n+1)*u*ma.pi)/(2*N))*ma.cos(((2*m+1)*v*ma.pi)/(2*N))

    return round((2/N)*x*y*fcos)



# Utile simplement pour voir le résultat, autrement j'utilise CDT() dans Cq(), pas CDTmatrice.
def CDTmatrice(matrice):
    l1 = []
    for i in range(len(matrice)):
        l2 = []
        for j in range(len(matrice)):
            l2.append(CDT(i,j, matrice))
        l1.append(l2)
        
    return l1
    


def Cq(matrice):
    l1 = []
    for i in range(len(matrice)):
        l2 = []
        for j in range(len(matrice)):
            l2.append(round(CDT(i,j, matrice)/PasQuantif[i][j]))
        l1.append(l2)
        
    return l1




def ZigZag(matrice):
    """
    Le principe est de lire en zig zag à partir des coordonées (ici sous forme de positions dans la liste de liste)

    Voilà à quoi ressemble la lecture en zig zag acutelle (1ere partie) :

    (0,0) 
    (1,0) (0,1) => Tuple pair
    (0,2) (1,1) (2,0) => Tuple impair
    (3,0) (2,1) (1,2) (0,3) => Tuple pair
    (0,4) (1,3) (2,2) (3,1) (4,0) => ...
    (5,0) (4,1) (3,2) (2,3) (1,4) (0,5) 
    (0,6) (1,5) (2,4) (3,3) (4,2) (5,1) (6,0)
    (7,0) (6,1) (5,2) (4,3) (3,4) (2,5) (1,6) (0,7)

    """
    coord = (0,0)
    liste = []
    max = 0
    for i in range(len(matrice)):
        max +=1
        for j in range(max):
            #Pour les "Tuple pair"
            if max%2 == 0 or max == 1:
                coord = (j,max-j-1)
            #Pour les "Tuple impair"
            else:
                coord = (i-j,j)
            liste.append(cq[coord[0]][coord[1]])
    
    """
    2eme boucle, pour l'autre moitié du carré :

    (1,7) (2,6) (3,5) (4,4) (5,3) (6,2) (7,1)
    (7,2) (6,3) (5,4) (4,5) (3,6) (2,7)
    (3,7) (4,6) (5,5) (6,4) (7,3)
    (7,4) (6,5) (5,6) (4,7)
    (5,7) (6,6) (7,5)
    (7,6) (6,7)
    (7,7)

    """    
    max = len(matrice)
    for i in range(len(matrice)):
        max -=1
        for j in range(max,0,-1):
            if max%2 == 0 or max == 1:
                coord = (j+i,len(matrice)-j)
            else:
                coord = (len(matrice)-j,j+i)
            liste.append(cq[coord[0]][coord[1]])

    return liste
            
            

def SansZigZag(matrice):
    """Fonction qui prend en valeur d'entrée une matrice a 2 dimension (liste double) et renvoie une matrice à 1 dimension (liste simple)

    Args:
        matrice (list): Matrice double

    Returns:
        liste: La liste simple (matrice 1 dimension)
    """
    liste = []
    for i in range(len(matrice)):
        for j in range(len(matrice[i])):
            liste.append(matrice[i][j])
    
    return liste



def RLC(liste):
    retliste = []
    j = 0
    while j < len(liste):
        char = liste[j]
        nb = 0
        while j != len(liste) and liste[j] == char:
            nb+=1
            j+=1
        print(nb, end=", ")
        retliste.append([nb,char])
    
    return retliste


def proba(liste):
    l = liste
    tot = 0
    prob = 0
    for i in liste:
        tot += i[0]
    
    for j in range(len(liste)):
        prob = round(liste[j][0]/tot, 3)
        l[j].append(prob)
    
    return l


def huffmann():
    pass





cq = Cq(Mbase)





go = True       
if go:            
    print("\n-------CDT DE LA MATRICE-------")
    print("Matrice obtenue : ")
    print(CDTmatrice(Mbase),"\n")
    print("Matrice voulu : ")
    print(ResVoulu,"\n")
    print("\n")


    print("-------CQ(u,v)-------")
    print(Cq(Mbase),"\n")
    print("\n")


    print("-------LECTURE ZIG ZAG-------")
    print(ZigZag(cq),"\n")
    print("Lecture sans Zig Zag")
    print(SansZigZag(cq))
    print("\n")


    print("-------RLC ET PROBA DU ZIG ZAG-------")
    print(proba(RLC(ZigZag(Mbase))),"\n")
    print("\n")


    print("-------RLC ET PROBA SANS ZIG ZAG-------") #Car j'ai l'impression de me tromper quelque part : Sans le ZigZag, le RLC à l'air plus efficace.
    print(proba(RLC(SansZigZag(cq))),"\n")