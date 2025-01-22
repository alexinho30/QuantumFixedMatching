from qiskit import QuantumCircuit
import matplotlib.pyplot as plt  
from math import log2, floor
from qiskit.quantum_info import random_statevector, Statevector
import utility as u

import random as r

def ROT(bitSequence, k, c):
    n = len(bitSequence)
    if(n & (n - 1) or k & (k -1)):
        print("errore bitSequence must have a length power of two\n")
        return

    print("bitSequence :" + str(bitSequence))

    qc = QuantumCircuit(n + 1)

    for i in range(n) : 
        if(bitSequence[i] == "1"): qc.x(i+1) 

    if(c): qc.x(0)

    for i in range(1, int(log2(n) - int(log2(k))+ 1)):
        for j in range(0, int(n/(k*2**i))):
            for x in range(k*j*2**i, k*j*2**i + k):
                qc.cswap(0, 1+x, 1+x + k*2**(i-1))


    qc.measure_all()

    results = u.run(qc, shots = 100)
    u.counts(results)


def M(x, y):
    n = len(x) 
    if(n != len(y)):
        print("x and y must have the same length")
        return 

    qc = QuantumCircuit(3*n, n)

    for i in range(n):
        if(x[i] == "1"): qc.x(i)
        if(y[i] == "1"): qc.x(i + n)
        qc.ccx(i, i + n, i + 2*n)
        qc.x(i)
        qc.x(i + n)
        qc.ccx(i, i+ n, i+ 2*n)

    qc.measure([i for i in range(2*n, 3*n)], list(range(n)))

    results = u.run(qc, shots = 100)
    u.counts(results)


def EXT(lamda, i):
    if(i < 0):
        print("error k must be greater  or equal than 0")
        return
    
    n = len(lamda)
    
    qc = QuantumCircuit(2*n - 2**i, n -2**i)

    for j in range(n):
        if(lamda[j] == "1"): qc.x(j)
    
    qc.barrier()
    
    for j in range(n - 2**i):
        qc.ccx(j, j + 2**i, j + n)

    qc.draw("mpl")

    qc.measure([j for j in range(n, 2*n -2**i)], [j for j in range(n - 2**i)])

    results = u.run(qc, shots = 100)
    u.counts(results)


def REVERSAL(x):
    n = len(x)
    if(n == 0):
        print("error you passed an empty string")
        return 

    qc = QuantumCircuit(n)

    for i in range(n):
        if(x[i] == "1"): qc.x(i)

    for i in range(0, floor(n/2)):
        qc.swap(i, n - i - 1)
    
    qc.measure_all()

    results = u.run(qc, shots = 100)
    u.counts(results)


def C(x, c):
    n = len(x)
    if(n == 0):
        print("error you passed an empty string")
        return 
    
    qc = QuantumCircuit(2*n + 1, n)

    if(not c): qc.x(0)

    for i in range(0,n):
        if(x[i] == "1"):
            qc.x(i + 1)
    
    for i in range(1, n + 1):
        qc.ccx(0, i, n + i)
    
    qc.measure([i for i in range(n + 1, 2*n + 1)], list(range(n)))
    results = u.run(qc, shots = 100)
    u.counts(results)

def BCJ(x, y, c):
    n = len(x)
    if(n != len(y)):
        print("error x and y must have the same length")
        return 
    
    qc = QuantumCircuit(3*n + 1, n)

    if(c): qc.x(0)

    for i in range(n):
        if(x[i] == "1"): qc.x(i + 1)
        if(y[i] == "1"): qc.x(i + n + 1)
    
    for i in range(1, n + 1):
        qc.mcx([0, i, i+n], i+ 2*n)
    
    qc.draw("mpl")
    plt.show()

    qc.measure([i for i in range(2*n + 1, 3*n + 1)], list(range(n)))
    results = u.run(qc, shots = 100)
    u.counts(results)
    
    
def DJ(x):
    n = len(x)
    if(n == 0):
        print("error x is an empty string")
        return
    
    qc = QuantumCircuit(n + 1, 1)

    for i in range(n):
        if(x[i] == "1"): qc.x(i)
    
    for i in range(n):
        qc.x(i)
    
    qc.mcx([i for i in range(n)], n)

    for i in range(n +1):
        qc.x(i)
    
    qc.measure(n, 0)
    results = u.run(qc, shots = 100)
    u.counts(results)
    







def lognImplementationFanOut(n, v):
    qc = QuantumCircuit(n)

    if(v): qc.x(0)

    i = 1
    while(i < n):
        for j in range(i):
            if(j + i < n): qc.cx(j, j+i)
        i= 2*i
    
    qc.measure_all()

    results = u.run(qc, 100)
   
    u.counts(results)


def mcxLognImplementation(n, v):
    qc = QuantumCircuit(2*n - 1)

    for i in range(n):
        if(v) : qc.x(i)

    k = n
    h = 0
    for i in range(int(log2(n))):
        for j in range(0, n, 2):
            qc.ccx(j + h, j+ h + 1, k)
            k+=1
        h+=n
        n = int(n/2)
        qc.barrier()
        
    qc.draw("mpl")
    plt.show()

    qc.measure_all()
    results = u.run(qc, 100)
    u.counts(results)




