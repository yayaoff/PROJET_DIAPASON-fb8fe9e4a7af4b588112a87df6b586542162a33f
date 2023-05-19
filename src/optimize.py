import subprocess
import time
import scipy.optimize
from math import sqrt

# Function optimize

r1 = 6e-3
r2 = 15e-3
e = 38e-3
l = 73e-3
f_target_0 = 1244.51
# f_target_1 = 622.25
output_file = "out_py.txt"

subprocess.run(['./opti_py', '1',output_file,str(r1),str(r2),str(e),str(l),'0'])
with open(output_file) as f:
    initial_freq_0 = float(f.readline())
print("Initial Freq 0 = " + str(initial_freq_0))

def fun_to_min_new(param_to_opt):
    if param_to_opt[2] >= r2:
        param_to_opt[2] -= r2
    # if param_to_opt[1] < param_to_opt[0]: 
        # param_to_opt[1] += param_to_opt[0]
        
    print("Computing Freq for r1={}, r2={}, e ={}, l={}".format(param_to_opt[2],r2,param_to_opt[1],param_to_opt[0]))
    subprocess.run( ['./opti_py', '1',output_file,str(param_to_opt[2]),str(r2),str(param_to_opt[1]),str(param_to_opt[0]),'0'])
    with open(output_file) as f:
        freq = float(f.readline())
    print('f = '+str(freq))
    return abs(f_target_0 - freq)

out = scipy.optimize.minimize(fun_to_min_new,(52e-3,12e-4,61e-4),bounds=[(1e-8,1e-1)],method='Nelder-Mead',tol=3,options = {'maxiter': 200000}).x
subprocess.run( ['./opti_py', '1',output_file,str(out[2]),str(r2),str(out[1]),str(out[0]),'1'])