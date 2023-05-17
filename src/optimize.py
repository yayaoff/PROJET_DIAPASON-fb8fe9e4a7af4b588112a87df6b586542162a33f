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
f_target_1 = 622.255
output_file = "out_py.txt"

subprocess.run(['./opti_for_py', '1',output_file,str(r1),str(r2),str(e),str(l)])
with open(output_file) as f:
    initial_freq_0 = float(f.readline())
    # initial_freq_1 = float(f.readline())
print("Initial Freq 0 = " + str(initial_freq_0))
# print("Initial Freq 1 = " + str(initial_freq_1))

def fun_to_min_1(param_to_opt):
    
    print("Computing Freq for r1={}, r2={}, e ={}, l={}".format(r1,param_to_opt[2],param_to_opt[1],param_to_opt[0]))
    subprocess.run( ['./opti_for_py', '1',output_file,str(r1),str(param_to_opt[2]),str(param_to_opt[1]),str(param_to_opt[0])])
    with open(output_file) as f:
        freq = float(f.readline())
    print('f = '+str(freq))
    return abs(f_target_0 - freq)

# scipy.optimize.minimize(fun_to_min_1,(1e-2,61e-3,58e-3),bounds=[(1e-7,4e-1)],method='Nelder-Mead')

def fun_to_min_new(param_to_opt):
    if param_to_opt[2] >= r2:
        param_to_opt[2] -= r2
    if param_to_opt[1] < param_to_opt[0]: 
        param_to_opt[1] += param_to_opt[0]
        
    print("Computing Freq for r1={}, r2={}, e ={}, l={}".format(param_to_opt[2],r2,param_to_opt[1],param_to_opt[0]))
    subprocess.run( ['./opti_for_py', '1',output_file,str(param_to_opt[2]),str(r2),str(param_to_opt[1]),str(param_to_opt[0])])
    with open(output_file) as f:
        freq = float(f.readline())
    print('f = '+str(freq))
    return abs(f_target_0 - freq)

scipy.optimize.minimize(fun_to_min_new,(52e-3,12e-4,61e-4),bounds=[(1e-8,1e-1)],method='Nelder-Mead')

def fun_to_min_2(param_to_opt):
    
    print("Computing Freq for r1={}, r2={}, e ={}, l={}".format(r1,param_to_opt[2],param_to_opt[1],param_to_opt[0]))
    subprocess.run( ['./opti_for_py', '2',output_file,str(r1),str(param_to_opt[2]),str(param_to_opt[1]),str(param_to_opt[0])])
    with open(output_file) as f:
        freq_0 = float(f.readline())
        freq_1 = float(f.readline())
    print('f_0 = '+str(freq_0))
    print('f_1 = '+str(freq_1))
    
    return sqrt(abs(f_target_0 - freq_1)**2 + abs(f_target_1- freq_0)**2)

# scipy.optimize.minimize(fun_to_min_2,(1e-3,61e-3,58e-2),bounds=[(1e-9,6e-1)],method='Nelder-Mead')
# scipy.optimize.minimize(fun_to_min_2,(1e-3,61e-3,58e-2),bounds=[(1e-9,6e-1)],method='Nelder-Mead',tol=1e-9)
# scipy.optimize.minimize(fun_to_min_2,(1e-3,61e-4,583e-3),bounds=[(1e-9,6e-1)],method='Nelder-Mead')