These two files are the code for the article "Digital Quantum Simulation and Circuit Learning for the Generation of Coherent States", which will be published on Entropy soon.
The environment required is as follow:
Python --3.6.5
qiskit --0.34.2
scipy --1.5.4
numpy --1.19.3
qutip --4.7.0

coherent_preparation_displacement_circuit.py : 
Preparation of coherent states by using DQS;
Corresponding to Chapter 2 of the article,
The main method in the file to prepare Coherent states is coherent_state(n,trotter_step,alpha)，
where n is the number of qubits (int), trotter_step is the number of trotter steps (int) and alpha is the displacement value (complex).
You can use this method to prepare any Coherent states with high fidelity by setting the appropriate trotter steps number and qubit number.

qiskit_all_scheme.py : 
Preparation of coherent states using three VQA schemes;
Here we have written all three Schemes used in the article in the file
You can reproduce our experimental data by changing the value of "scheme"
You can also change the circuit to build your own Ansatz

Please do not hesitate to contact us by email when having any questions.

Email: L_ruilin@shu.edu.cn
