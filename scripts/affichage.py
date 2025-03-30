### =================================
### Import des librairies nécessaires
### =================================
from matplotlib import pyplot as plt
from matplotlib import gridspec as gridspec
import numpy as np
from scipy.special import assoc_laguerre, genlaguerre, factorial
plt.style.use("article.mplstyle")


### ================================================
### Lecture du fichier pour la fonction d'onde radiale
### ================================================
Nom_fichier = "data/valeur_radiale.out"
fichier = open(Nom_fichier, "r")
Lignes = fichier.readlines()

# Extraction des paramètres physiques
E, Z, l, n, dr = float(Lignes[0].split(";")[0][2:]), float(Lignes[0].split(";")[1][2:]), \
                 float(Lignes[0].split(";")[2][2:]), float(Lignes[0].split(";")[3][2:]), \
                 float(Lignes[0].split(";")[4][3:])

# Dictionnaire pour la notation spectroscopique
Dico_l = {0:"s", 1:"p", 2:"d", 3:"f", 4:"g"}

# Lecture des données numériques (r, u(r))
X, Y = [], []
for i in range(len(Lignes)):
    if i >= 3:
        ligne = Lignes[i].split(";")
        X.append(float(ligne[1]))  # Valeurs de r
        Y.append(float(ligne[0]))  # Valeurs de u(r)

# Conversion en tableaux numpy et calcul de la densité
A_r = np.array(X)  # Tableau des rayons
A_u = np.array(Y)  # Tableau des valeurs de la fonction d'onde
A_dense = A_u**2/A_r**2  # Densité de probabilité radiale
norm = np.cumsum(A_dense)[-1]  # Normalisation


### ================================================
### Calcul théorique de la densité de probabilité
### ================================================
def densite_probabilite(r, n, l, Z):
    """
    Calcule la fonction d'onde radiale théorique pour un hydrogénoïde

    Paramètres:
    r (array): Points radiaux
    n (int): Nombre quantique principal
    l (int): Nombre quantique orbital
    Z (float): Charge nucléaire

    Retourne:
    array: Fonction d'onde radiale normalisée
    """
    a0 = 1  # Rayon de Bohr en unités atomiques
    r = r / a0

    # Calcul du facteur de normalisation et de la partie radiale
    rho = np.sqrt((2 * Z / (n * a0))**3 * factorial(n - l - 1) / (2 * n * factorial(n + l)**3))
    rho *= np.exp(- Z * r / (n * a0)) * (2 * Z * r / (n * a0))**(l)

    # Calcul du polynôme de Laguerre généralisé
    L = genlaguerre(n - l - 1, 2 * l + 1)
    rho *= L(2 * Z * r / (n * a0))

    return rho

### ==============================================
### Calcul des densités théoriques et simulation
### ==============================================
# Calcul de la densité théorique et sa normalisation
A_dense_theo = densite_probabilite(A_r,n,l,Z)**2
norm_theo = np.cumsum(A_dense_theo)[-1]


### ==============================================
### Affichage des courbes de densité radiale
### ==============================================
plt.figure("équation_radiale")
plt.ticklabel_format(axis='y', style="scientific", scilimits=(0, 0), useMathText=True)
# Configuration du titre et des axes avec unités atomiques
plt.title("Z="+str(Z)+ ", E="+ str(E) + " [E$_h$]" +", "+str(int(n))+Dico_l[l])
plt.xlabel("Rayon [a$_0$]")
plt.ylabel("$|R(r)|^2 [1]$")

# Tracé des courbes normalisées
plt.plot(A_r,A_dense/norm,label="Simulation")
plt.plot(A_r,A_dense_theo/norm_theo,label="Théorie",color="blue",linestyle="dashed", dashes=(10,20))
plt.legend()
print("u(0)^2=",A_dense[-1]/norm)


### ==============================================
### Affichage de la différence théorie-simulation
### ==============================================
plt.figure("différence")
plt.ticklabel_format(axis='y', style="scientific", scilimits=(0, 0), useMathText=True)
plt.title("Z="+str(Z)+ ", E="+ str(E) + " [E$_h$]" +", "+str(int(n))+Dico_l[l])
plt.xlabel("Rayon [a$_0$]")
plt.ylabel("$|R_{simu}(r)|^2 - |R_{théo}(r)|^2 [1]$")
plt.plot(A_r,A_dense/norm-A_dense_theo/norm_theo,label="Simulation - Théorie",color="green")
plt.legend()


### ==============================================
### Lecture des données pour la minimisation
### ==============================================
Nom_fichier = "data/valeur_rfixe.out"
fichier = open(Nom_fichier,"r")
Lignes = fichier.readlines()

# Extraction des paramètres physiques pour la minimisation
Z,l,dr = float(Lignes[0].split(";")[0][2:]), \
         float(Lignes[0].split(";")[1][2:]), \
         float(Lignes[0].split(";")[2][3:])

# Configuration des niveaux d'énergie à afficher
n_min,n_max = int(l+1),6

# Lecture des données (E, S(E))
X,Y = [],[]
for i in range(len(Lignes)):
    if i >= 3:
        ligne = Lignes[i].split(";")
        X.append(float(ligne[1]))  # Valeurs de E
        Y.append(float(ligne[0]))  # Valeurs de S(E)
A_E = np.array(X)
A_S = np.array(Y)


### ==============================================
### Fonctions d'affichage des résultats
### ==============================================
def deux_plot():
    """
    Affiche les résultats sur deux graphiques côte à côte
    pour une meilleure visualisation des différentes zones d'énergie
    """
    fig = plt.figure(figsize=(8,6),num="équation_r_fixe")
    gs = gridspec.GridSpec(1,2)
    ax1,ax2 = fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1])
    fig.suptitle("$Z="+str(Z)+", l="+str(l)+"$",fontsize=25)

    ax1.set_xlabel("Energie [E$_h$]")
    ax1.set_ylabel(r"$|U(0,E)|^2dr$ [1]",fontsize=20)

    ax2.set_xlabel("Energie [E$_h$]")
    #ax2.set_ylabel(r"$\frac{|U(0)|^2dr}{max(|U(r)|^2dr)}$ [1]")

    split_index = int((np.abs(A_E - -Z**2/(2*3**2)).argmin() + np.abs(A_E - -Z**2/(2*4**2)).argmin())/2)+100
    ax1.set_xlim((np.min(A_E[:split_index]),np.max(A_E[:split_index])))
    ax2.set_xlim((np.min(A_E[split_index:]),np.max(A_E[split_index:])))

    ax1.set_ylim(0,1.1)
    ax2.set_ylim(0,1.1)

    ax1.plot(A_E[:split_index],A_S[:split_index])
    ax2.plot(A_E[split_index:],A_S[split_index:])
    max_dense = np.max(A_S)
    for n in range(n_min,4):
        if n>l:
            if n==l+1:
                ax1.plot([-Z*Z/(2*n**2),-Z*Z/(2*n**2)],[0,max_dense],"--",color = 'blue',label=r"E$_n$ = $\frac{-Z^2}{2n^2}$")
                ax1.text(-Z*Z/(2*n**2),max_dense+0.025,str(n),fontsize=16,zorder=10,ha="center")
            else:
                ax1.plot([-Z*Z/(2*n**2),-Z*Z/(2*n**2)],[0,max_dense],"--",color = 'blue')
                ax1.text(-Z*Z/(2*n**2),max_dense+0.025,str(n),fontsize=16,zorder=10,ha="center")

    for n in range(4,n_max+1):
        if n>l:
            if n==l+1:
                ax2.plot([-Z*Z/(2*n**2),-Z*Z/(2*n**2)],[0,max_dense],"--",color = 'blue',label=r"E$_n$ = $\frac{-Z^2}{2n^2}$")
                ax2.text(-Z*Z/(2*n**2),max_dense+0.025,str(n),fontsize=16,zorder=10,ha="center")
            else:
                ax2.plot([-Z*Z/(2*n**2),-Z*Z/(2*n**2)],[0,max_dense],"--",color = 'blue')
                ax2.text(-Z*Z/(2*n**2),max_dense+0.025,str(n),fontsize=16,zorder=10,ha="center")

    plt.ticklabel_format(axis='x', style="scientific", scilimits=(0, 0), useMathText=True)
    ax1.ticklabel_format(axis='x', style="scientific", scilimits=(0, 0), useMathText=True)
    ax2.ticklabel_format(axis='x', style="scientific", scilimits=(0, 0), useMathText=True)
    ax1.legend(fontsize=18)
    plt.show()

def un_plot():
    """
    Affiche les résultats sur un seul graphique
    """
    fig = plt.figure(figsize=(8,6),num="équation_r_fixe")
    gs = gridspec.GridSpec(1,1)
    ax1 = fig.add_subplot(gs[0,0])
    ax1.set_title("$Z="+str(Z)+", l="+str(l)+'$')

    ax1.set_xlabel("Energie [E$_h$]")
    ax1.set_ylabel(r"$|U(0,E)|^2dr}$ [1]",fontsize=20)

    ax1.set_xlim((np.min(A_E),np.max(A_E)))

    ax1.set_ylim(0,1.1)

    ax1.plot(A_E,A_S)
    max_dense = np.max(A_S)
    for n in range(n_min,n_max+1):
        if n>l:
            if n==l+1:
                ax1.plot([-Z*Z/(2*n**2),-Z*Z/(2*n**2)],[0,max_dense],"--",color = 'blue',label=r"E$_n$ = $\frac{-Z^2}{2n^2}$")
                ax1.text(-Z*Z/(2*n**2),max_dense+0.025,str(n),fontsize=16,zorder=10,ha="center")
            else:
                ax1.plot([-Z*Z/(2*n**2),-Z*Z/(2*n**2)],[0,max_dense],"--",color = 'blue')
                ax1.text(-Z*Z/(2*n**2),max_dense+0.025,str(n),fontsize=16,zorder=10,ha="center")

    plt.ticklabel_format(axis='x', style="scientific", scilimits=(0, 0), useMathText=True)
    ax1.legend(fontsize=18)
    plt.show()

# Choix de la fonction d'affichage
deux_plot()
#un_plot()