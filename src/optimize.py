import subprocess
import time
import scipy.optimize

# Function optimize

r1 = 6e-3
r2 = 11e-3
e = 38e-3
l = 82e-3
f_target = 1244.51
output_file = "out_py.txt"

subprocess.run(['./opti_for_py', '1',output_file,str(r1),str(r2),str(e),str(l)])
with open(output_file) as f:
    initial_freq = float(f.readline())
print("Initial Freq = " + str(initial_freq))

def fun_to_min(param_to_opt):
    
    print("Computing Freq for r1={}, r2={}, e ={}, l={}".format(r1,param_to_opt[2],param_to_opt[1],param_to_opt[0]))
    subprocess.run( ['./opti_for_py', '1',output_file,str(r1),str(param_to_opt[2]),str(param_to_opt[1]),str(param_to_opt[0])])
    with open(output_file) as f:
        freq = float(f.readline())
    print('f = '+str(freq))
    return abs(f_target - freq)

scipy.optimize.minimize(fun_to_min,(1e-2,61e-3,58e-3),bounds=[(1e-7,4e-1)],method='Nelder-Mead')