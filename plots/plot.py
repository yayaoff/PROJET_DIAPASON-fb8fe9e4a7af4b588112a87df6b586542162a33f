import numpy as np
import matplotlib.pyplot as plt

# Load the matrix from the CSV file
# matrix=np.genfromtxt('A.csv', delimiter=' ')
# matrix_band = np.genfromtxt('A_band.csv', delimiter=' ')
time_creuse= np.genfromtxt('matrix.csv', delimiter=',', names=True, dtype=float)
time_band = np.genfromtxt('matrix_band.csv', delimiter=',', names=True,dtype=float)
time_creuse_half= np.genfromtxt('matrix_half.csv', delimiter=',', names=True, dtype=float)
time_band_half = np.genfromtxt('matrix_band_half.csv', delimiter=',', names=True,dtype=float)

def plot_mask(matrix,title):
    # Create a binary mask of the non-zero elements in the matrix
    mask = matrix != 0

    # Plot the mask as a black-and-white image
    plt.imshow(mask, cmap='binary')
    plt.title(title)
    #pour sauvegarder les images, utiliser plt.savefig
    plt.show()

def plot_time(data1, data2, data3, data4,title):
    size=[]
    creuse=[]
    bande=[]
    creuse_h=[]
    bande_h=[]

    
    for i in range(len(data1)):
        size.append(data1[i][0])
        creuse.append(data1[i][1])
        bande.append(data2[i][1])
        creuse_h.append(data3[i][1])
        bande_h.append(data4[i][1])

    n3=np.zeros(len(size))
    n_bande=np.zeros(len(size))
    for i in range(len(n3)):
        n3[i]=size[i]**3+size[i]**2
        n_bande[i]=size[i]**2

    #pour sauvegarder les images, utiliser plt.savefig
    #plot en normal
    #---------------
    # plt.plot(size, creuse, label='Matrice Creuse')
    # plt.plot(size, bande, label='Matrice Bande')
    # plt.plot(size, creuse_h, label='Demi Diapason Creuse')
    # plt.plot(size, bande_h, label='Demi Diapason, Bande')
    # plt.xlabel('Taille de Matrice[n]')
    # plt.ylabel('Temps de résolution [s]')

    #plot en echelle logarithmique
    #-----------------------------
    # plt.loglog(size, creuse/np.linalg.norm(creuse), label='Matrice Creuse')
    # plt.loglog(size, creuse_h/np.linalg.norm(creuse_h), label='Demi Diapason Creuse')
    # plt.loglog(size, n3/np.linalg.norm(n3),label='Courbe Théorique Creuse')
    plt.loglog(size, bande_h/np.linalg.norm(bande_h), label='Demi Diapason, Bande')
    plt.loglog(size, bande/np.linalg.norm(bande), label='Matrice Bande')
    plt.loglog(size, n_bande/np.linalg.norm(n_bande),label='Courbe Théorique Bande')
    plt.xlabel('Taille de Matrice[log(n)]')
    plt.ylabel('Temps de résolution [log(s)]')
    plt.title(title)
    plt.legend()
    plt.show()

#commenter/decommenter en fction de ce qu'on veut
# plot_mask(matrix, "Masque de la matrice creuse")
# plot_mask(matrix_band,"Masque de la matrice bande")
plot_time(time_creuse,time_band,time_creuse_half,time_band_half,"Comparaison du temps d'execution avec la courbe théorique")