import subprocess
import time

target1 = 770
target2 = 1540

r1 = 6e-3
r2 = 13e-3
e = 38e-3
l1 = 0.046240760
borne_infl2 = 0.05
borne_supl2 = 0.1
borne_infl1 = 0.04
borne_supl1 = 0.06

ok = False
tol = 0.1
nb_iter = 0

###premier paramètre
t1 = time.time()

while(ok == False):
    nb_iter+=1
    new_bornel2 = (borne_infl2+borne_supl2)/2
    subprocess.run(['./opti_for_py', '1', 'out.txt',str(r1),str(r2),str(e),str(l1),str(new_bornel2)])
    with open('out.txt', 'r') as f:
        ligne = f.readline().strip()
        freqs = [float(num) for num in ligne.split()]
        freq = freqs[0]
    #print('\n \n Fréquence = ', freq)
    if(abs(freq-target1)<tol):
        print("\n \n longueure optimale trouvée ! ", new_bornel2)
        ok = True
    if(freq-target1 > 0):
        borne_infl2 = new_bornel2
    else:
        borne_supl2 = new_bornel2
    print("\n \n \n")
ok = False
while(ok == False):
    nb_iter+=1
    new_bornel1 = (borne_infl1+borne_supl1)/2
    subprocess.run(['./opti', '2', 'out.txt',str(r1),str(r2),str(e),str(new_bornel1),str(new_bornel2)])
    with open('out.txt', 'r') as f:
        ligne = f.readline().strip()
        freqs = [float(num) for num in ligne.split()]
        freq = freqs[1]
    #print('\n \n Fréquence = ', freq)
    if(abs(freq-target2)<tol):
        print("\n \n longueure optimale 2 trouvée ! ", new_bornel1)
        ok = True
    if(freq-target2 > 0):
        borne_infl1 = new_bornel1
    else:
        borne_supl1 = new_bornel1
    print("\n \n \n")
t2 = time.time()
temps = t2-t1
minutes, seconds = divmod(temps, 60) # convert timestamp to minutes and seconds
time_str = time.strftime('Temps de l\'optimisation :%M min %S sec ', time.gmtime(temps)) # format minutes and seconds as MM:SS

subprocess.run(['./opti', '2', 'out.txt',str(r1),str(r2),str(e),str(new_bornel1),str(new_bornel2)])
with open('out.txt', 'r') as f:
    ligne = f.readline().strip()
    freqs = [float(num) for num in ligne.split()]
print("fin de la bissection ! \n ")
print(time_str) 
print("\n")
print("longueur extérieure: ", new_bornel2)
print("\n")
print("longueur intérieure: ", new_bornel1)
print("\n")
print("première fréquence: ",freqs[0])
print("\n")
print("deuxième fréquence: ",freqs[1])
print("\n")
print("Nombre d'itérations",nb_iter)