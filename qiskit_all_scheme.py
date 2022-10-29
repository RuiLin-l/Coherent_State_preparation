# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 19:53:18 2021

@author: Administrator
"""
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, transpile,execute,Aer
from qiskit.quantum_info import Statevector,state_fidelity
from qiskit.result import Result
from qiskit import *
import numpy as np
import math
from qiskit.quantum_info import Statevector,state_fidelity
from qutip import *


def VQA_CHOERENT_SCHEME_A(n,parameters,depth):
    '''
    Circuit of Scheme a;
    n:qubit numbers;
    depth:circuit depth;
    '''
    circuit = QuantumCircuit(n, n) 


    for K in range(depth):
        for i in range(n):
        
            circuit.rx(parameters[3*i+4*K*n],i)
            circuit.rz(parameters[3*i+4*K*n+1],i)
            circuit.rx(parameters[3*i+4*K*n+2],i)
          
        circuit.cry(parameters[4*K*n+3*n], 0, 1, label=None, ctrl_state=None)
        circuit.cry(parameters[4*K*n+3*n+1], 1, 2, label=None, ctrl_state=None)
        circuit.cry(parameters[4*K*n+3*n+2], 2, 3, label=None, ctrl_state=None)
        circuit.cry(parameters[4*K*n+3*n+3], 3, 0, label=None, ctrl_state=None)

 
    circuit.save_statevector() 
    circuit.measure_all() 
    simulator = Aer.get_backend("aer_simulator")
    circuit = transpile(circuit,simulator)
    result = simulator.run(circuit).result()
    statevector_reverse = Statevector(result.get_statevector(circuit)).reverse_qargs()
    result = statevector_reverse
    return result

def VQA_CHOERENT_SCHEME_B(n,parameters,depth):
    '''
    Circuit of Scheme b;
    n:qubit numbers;
    depth:circuit depth;
    '''
    circuit = QuantumCircuit(n, n) 
    parameters_1 = parameters[:3*n]
    parameters_2 = parameters[3*n:]
    for i in range(n):
        circuit.rx(parameters_1[3*i],i)
        circuit.rz(parameters_1[3*i+1],i)
        circuit.rx(parameters_1[3*i+2],i)

    for k in range(depth):
        
        for i in range(1,n):
            circuit.cx(i-1,i)     
        for j in range(n):

            circuit.ry(parameters_2[j+k*n],j)
        

        
    circuit.save_statevector() 
    circuit.measure_all() 
    simulator = Aer.get_backend("aer_simulator")
    circuit = transpile(circuit,simulator)
    result = simulator.run(circuit).result()
    statevector_reverse = Statevector(result.get_statevector(circuit)).reverse_qargs()
    result = statevector_reverse
    return result
def VQA_CHOERENT_SCHEME_C(n,parameters_1,depth):
    '''
    Circuit of Scheme c;
    n:qubit numbers;
    depth:circuit depth;
    '''
    circuit = QuantumCircuit(n, n) 
    for i in range(depth):

        circuit.rx(parameters_1[15*i],0)
        circuit.rx(parameters_1[15*i+1],1)
        circuit.cx(0,1)
        circuit.rz(parameters_1[15*i+2],1)
        circuit.cx(0,1)
        circuit.rz(parameters_1[15*i+3],0)
        circuit.rz(parameters_1[15*i+4],1)        
        circuit.rx(parameters_1[15*i+5],1)
        circuit.rx(parameters_1[15*i+6],2)
        circuit.cx(1,2)
        circuit.rz(parameters_1[15*i+7],2)
        circuit.cx(1,2)
        circuit.rz(parameters_1[15*i+8],1)
        circuit.rz(parameters_1[15*i+9],2)
     
        circuit.rx(parameters_1[15*i+10],2)
        circuit.rx(parameters_1[15*i+11],3)
        circuit.cx(2,3)
        circuit.rz(parameters_1[15*i+12],3)
        circuit.cx(2,3)
        circuit.rz(parameters_1[15*i+13],2)
        circuit.rz(parameters_1[15*i+14],3)   

    circuit.save_statevector() 
    circuit.measure_all() 
    simulator = Aer.get_backend("aer_simulator")
    circuit = transpile(circuit,simulator)

    result = simulator.run(circuit).result()
    statevector_reverse = Statevector(result.get_statevector(circuit)).reverse_qargs()
    result = statevector_reverse

    return result


if __name__ == '__main__':
    
    exact_value = qutip.coherent(16,1+1j,method = 'analytic')


    def get_statefidelity(a):        

        fidelity = state_fidelity(exact_value[0:],a,validate = False)
        print("fidelity",fidelity)
        infidelity = 1-fidelity
        
        return infidelity
    
    def cost_func(parameters):
        if scheme == "SCHEME_A":
            
            a= VQA_CHOERENT_SCHEME_A(n,parameters,depth)

        elif scheme == "SCHEME_B":
            a= VQA_CHOERENT_SCHEME_B(n,parameters,depth)

        elif scheme == "SCHEME_C":
            a= VQA_CHOERENT_SCHEME_C(n,parameters,depth)

        else:
            print("no scheme")
        fidelity = get_statefidelity(a)   
        return fidelity 
    n = 4
    depth = 6
    scheme = "SCHEME_B"

    
    parameters_A = [1]*(4*n*depth) #initial parameters of Scheme a
    parameters_B = [1]*(3*n+4*depth) #initial parameters of Scheme b
    parameters_C = [1]*(15*depth) #initial parameters of Scheme c
    


    if scheme == "SCHEME_A":
        parameters = parameters_A
    elif scheme == "SCHEME_B":
        parameters = parameters_B
    elif scheme == "SCHEME_C":
        parameters = parameters_C
    else:
        print("no scheme")

    res = minimize(cost_func,parameters,method='SLSQP',options= {'maxiter':5000 })
    print(res)


    