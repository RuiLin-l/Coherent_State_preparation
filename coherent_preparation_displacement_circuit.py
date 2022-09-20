
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, transpile,execute,Aer
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram,plot_state_city,plot_bloch_multivector,plot_state_qsphere
from qiskit.quantum_info import Statevector,state_fidelity
from qiskit.result import Result
from qiskit import *
import numpy as np
import math
from qiskit.quantum_info import Statevector,state_fidelity
from qutip import *

#eng=MainEngine(backend = Simulator(gate_fusion=True))

def aa(n):
    a = []

    for i in range(2**n):

        a_int = []
        for j in range(2**n):
            if j == i+1:
                a_int.append(np.sqrt(j))
            else:
                a_int.append(0)
        a.append(a_int)

    return a
def aadagger(n):
    a = []
    for i in range(2 ** n):

        a_int = []
        for j in range(2 ** n):
            if j == i - 1:
                a_int.append(np.sqrt(i))
            else:
                a_int.append(0)
        a.append(a_int)
    return a
def decompose_item(qubit_n,row,col):

    n = 2**qubit_n
    result_str=""
    while (n !=1):
        if row <=n/2:
            if col <=n/2:
                result_str+="1"
            elif col >n/2:
                result_str+="2"
        elif row >=n/2:
            if col <=n/2:
                result_str+="3"
            elif col >n/2:
                result_str+="4"
        n = n/2
        if row >n:
            row = row-n
        else:
            row = row
        if col >n:
            col = col-n
        else:
            col = col
    return result_str
def matrix_item(matrix):

    row=0
    matrix_item_dict={}
    for row_item in matrix:
        row+=1
        col=0
        for col_item in row_item:
            col+=1
            if col_item !=0:
                matrix_item_dict[f"{row},{col}"] = col_item
    return matrix_item_dict
def decompose(qubit_n,matrix):

    matrix_item_dict = matrix_item(matrix)
    matrix_item_dict=dict(sorted(matrix_item_dict.items(), key=lambda x: x[1], reverse=False))
    dict2={}
    for coordinate,val in matrix_item_dict.items():
        row = int(coordinate.split(",")[0])
        col = int(coordinate.split(",")[1])
        cof = val
        decompose_item_result = decompose_item(qubit_n,row,col)
        dict2[decompose_item_result] = val
    return dict2
def string2pauli(str1): 

    list_all = []
    for i in range(len((str1))):
        if str1[i] == "1":
            list_all.append(("I","Z"))
        elif str1[i] == "2":
            list_all.append(("X","-iY"))
        elif str1[i] == "3":
            list_all.append(("X","iY"))
        elif str1[i] == "4":
            list_all.append(("I","-Z"))
    return list_all
def merge(i,j):

    res = [] 
    for p in (i,j): 
        if isinstance(p,tuple): 
            res.extend(p) 
        else: res.append(p) 
    return tuple(res)

def combineN(*args):

    target = args[0] 
    for li in args[1:]:
        tmp = [] 
        for i in target: 
            for j in li: 
                s = "".join((merge(i,j)))
                tmp.append(s)      
        target = tmp 
    return target
def change_style(str1):

    c = string2pauli(str1)
    sigma2pauli = combineN(*c)
    list_change_i=[] 
    for i in sigma2pauli:        
        if i.count("i")>1:
            if i.count("i")%2 == 0: 
                i=(int((i.count("i"))/2))*"-"+i
                i=i.replace("i","")                
            elif i.count("i")%2 != 0:
                i=(int((i.count("i"))/2))*"-"+i
                i=i.replace("i","")                
                i = "i"+i
        elif i.count("i")==1:
            i = i.replace("i","")
            i = "i"+i
        else:
            i = i
        list_change_i.append(i)
    list_change_symbol=[]
    for i in list_change_i:
        if i.count("-") == 1:
            i=i.replace("-","")
            i = "-"+i
        elif i.count("-")>1:
            if i.count("-")%2 == 0:
                i = i.replace("-","")
            else:
                i = i.replace("-","")
                i = "-"+i
        else:
            i=i 
        list_change_symbol.append(i)
    return list_change_symbol

def matrix2dict(n,matrix):

    dict1=decompose(n,matrix)
    dict3={}   
    for keys in dict1.keys():
        list_after_keys =change_style(keys)
        for item in list_after_keys:
            if item.replace("-","") in dict3.keys():
                if item[0] == "-":
                    item = item.replace("-","")
                    dict3[item] = dict3[item]-dict1[keys]
                else:
                    dict3[item] = dict3[item]+dict1[keys]
            else:
                if item[0] == "-": 
                    item = item.replace("-","")
                    dict3[item] =  -dict1[keys]
                else:
                    dict3[item] =  dict1[keys]
    dict4={}
    for keys in dict3.keys():
        if dict3[keys]  != 0.0:            
            dict4[keys] = dict3[keys]
    dict5 = {}
    
    for keys in dict4.keys():
        if dict4[keys].imag != 0 and keys[0] == "i":
            dict5[keys.replace("i","")] = dict4[keys].imag
        elif dict4[keys] >= 0 and dict4[keys] >=1e-6:
            dict5[keys] = dict4[keys]
        elif dict4[keys] <= 0 and dict4[keys] <=-1e-6:
            dict5[keys] = dict4[keys]
    return(dict5)
def key_circuit(circuit,dict1,coef):

    
    for str1 in dict1.keys():

        list_str1 = list(str1)
        for index,item in enumerate(list_str1):    
            if item =="X":
                circuit.h(index)
        for index,item in enumerate(list_str1):
            if item =="Y":
                circuit.rx((np.pi/2),index)                    
        str2 = str1.replace("X","Z").replace("Y","Z")
        list_str2=list(str2)
        index_last = str2.rfind("Z")
        for index,item in enumerate(list_str2):    
            if item =="Z" and index != index_last:
                circuit.cx(index,index_last)
        circuit.rz((dict1[str1]*coef),index_last) 
        for index,item in enumerate(list_str2):
            index_last = str2.rfind("Z")
            if item =="Z" and index != index_last:
                circuit.cx(index,index_last)
 
        for index,item in enumerate(list_str1):    
            if item =="Y":
                circuit.rx((-np.pi/2),index)
        for index,item in enumerate(list_str1):
            if item =="X":  
                circuit.h(index)
    return circuit


def coherent_state(n,trotter_step,alpha):
    
    real_part = alpha.real
    imag_part = alpha.imag
    a = aa(n)
    a_dagger = aadagger(n)
    c = np.dot(a,a)
    Z1 = 1j*np.array(np.array(a_dagger)-np.array(a))
    Z2 = -(np.array(np.array(a)+np.array(a_dagger)))
    dict1 = matrix2dict(n,Z1)
    dict2 = matrix2dict(n,Z2)

    real_part_item = real_part / trotter_step
    imag_part_item = imag_part / trotter_step
    circuit = QuantumCircuit(n, n)  
    coef1 = real_part_item/(2**(n-1))
    coef2 = imag_part_item/(2**(n-1))
    for j in np.arange(trotter_step):
        circuit = key_circuit(circuit,dict1,coef1)
        circuit = key_circuit(circuit,dict2,coef2)
    circuit.save_statevector() 
    circuit.measure_all() 
    simulator = Aer.get_backend("aer_simulator")
    circuit = transpile(circuit,simulator)
    result = simulator.run(circuit).result()
    statevector = result.get_statevector(circuit)
    statevector_reverse = Statevector(result.get_statevector(circuit)).reverse_qargs()
    result = statevector_reverse
    return result 


if __name__ == '__main__':

    n = 4
    print(coherent_state(n,10,1))

  
    
