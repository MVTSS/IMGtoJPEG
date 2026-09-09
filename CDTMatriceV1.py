"""
Réduction d'une matrice avec la transformée en cosinus discret
Matisse A.

v1
"""

from math import cos, pi, sqrt

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
        return 1/sqrt(2) 
    else:
        return 1

def CDT(u,v, matrice):
    fcos = 0
    N = len(matrice)
    for n in range(N):
        for m in range(N):
            #f(n,m, matrice) to (f(n,m, matrice) - 128) because the range of a pixel in JPEG is from [-128,127] instead of [0,255]
            fcos += ((f(n,m, matrice) - 128)*cos(((2*n+1)*u*pi)/(2*N))*cos(((2*m+1)*v*pi)/(2*N)))   

    return round((2/N)*c(u)*c(v)*fcos)



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
    
    Idea took from Zig Zag Matrix Diagonal Traversal online
    """
    N = len(matrice)
    resultat = []
    for s in range(2 * N - 1):
        if s % 2 == 0:
            # Remonter en diagonale
            r = min(s, N - 1)
            c = s - r
            while r >= 0 and c < N:
                resultat.append(matrice[r][c])
                r -= 1
                c += 1
        else:
            # Descendre en diagonale
            c = min(s, N - 1)
            r = s - c
            while c >= 0 and r < N:
                resultat.append(matrice[r][c])
                r += 1
                c -= 1
    return resultat
            
            

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







def print_matrix(title, matrice):
    print(f"\n=== {title} ===")
    for line in matrice:
        print(" ".join(f"{val:4d}" for val in line))
           
cq = Cq(Mbase)
cdt_obtenue = CDTmatrice(Mbase)
zigzag_cq = ZigZag(cq)
sans_zigzag_cq = SansZigZag(cq)

print_matrix("CDT OBTENUE", cdt_obtenue)
print_matrix("RÉSULTAT VOULU (TD)", ResVoulu)
print_matrix("MATRICE QUANTIFIÉE Cq", cq)

print("\n=== COMPARISON RLC (ZIGZAG vs SANS ZIGZAG) ===")
print(f"Nombre de paires RLC AVEC ZigZag : {len(RLC(zigzag_cq))}")
print(f"Nombre de paires RLC SANS ZigZag : {len(RLC(sans_zigzag_cq))}")

print("\n=== RLC + PROBABILITÉS (ZIGZAG) ===")
print("Format: [Occurrences, Valeur, Probabilité]")
for item in proba(RLC(zigzag_cq)):
    print(item)
