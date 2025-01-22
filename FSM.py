from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from math import log2, floor, ceil
import utility as u


def initializeCircuit(x, y, patternLength):
    binD = bin(patternLength)[2:]

    numIter = floor(log2(patternLength)) + 1
    n = len(x)
    m = len(binD)

    r = ClassicalRegister(1, "result")

    QuantumRegisters = []

    QuantumRegisters.append(QuantumRegister(n, name="x"))
    QuantumRegisters.append(QuantumRegister(n, "y"))
    QuantumRegisters.append(QuantumRegister(m, "d"))

    for i in range(numIter):
        QuantumRegisters.append( QuantumRegister(n - 2**i + 1, name="λ" + str(i)))

    for i in range(-1, numIter):
        QuantumRegisters.append(QuantumRegister(n, name="D" + str(i)))
    
    QuantumRegisters.append(QuantumRegister(1, "v"))

    
    qc = QuantumCircuit(*QuantumRegisters, r)


    print("numQbits used:")
    print(2*n + m + m*(n) + (m+1)*n + 1)

    #initializing input strings
    for i in range(n):
        if(x[i] == "1"): qc.x(i)
        if(y[i] == "1"): qc.x(n+i)

    #initializing binary version of d
    for i in range(len(binD)):
        if(binD[i] == "1"): qc.x(2*n + i)

    return qc 

def REV(dLen):
    qc = QuantumCircuit(dLen)

    for i in range(0, floor(dLen/2)):
        qc.swap(i, dLen - i - 1)
    
    qc = qc.to_gate(label="REV")
    return qc 

def M(n):
    qc = QuantumCircuit(3*n)

    for i in range(n):
        qc.ccx(i, i + n, i + 2*n)
        qc.x(i)
        qc.x(i + n)
        qc.ccx(i, i+ n, i+ 2*n)
        qc.x(i)
        qc.x(i + n)
    
    qc = qc.to_gate(label="M")
    return qc 

def EXT(n, i):

    qc = QuantumCircuit(2*n - 2**i)

    for j in range(n -2**i):
        qc.ccx(j, j + 2**i, j + n)

    qc = qc.to_gate(label="EXT" + str(i))

    return qc 

def ID(n):

    qc = QuantumCircuit(n)

    for i in range(n):
        qc.x(i)
    
    qc = qc.to_gate(label="ID")
    return qc 

def ROT(n, k):

    qc = QuantumCircuit(n +1)

    for i in range(1, int(log2(n) - int(log2(k))+ 1)):
        for j in range(0, int(n/(k*2**i))):
            for x in range(k*j*2**i, k*j*2**i + k):
                qc.cswap(0, x + 1, 1 + x + k*2**(i-1))
    
    qc = qc.to_gate(label="ROT" + str(k))
    return qc 

def BCJ(lamdaDim, n):

    qc = QuantumCircuit(lamdaDim + 2*n + 1)

    for i in range(lamdaDim):
        qc.mcx([0, i +1, lamdaDim + i + 1], lamdaDim + n + i)
    
    qc = qc.to_gate(label="BCJ")
    return qc 

def C(n):

    qc = QuantumCircuit(2*n + 1)

    qc.x(0)

    for i in range(n):
        qc.ccx(0, i + 1, i +n + 1)
    
    qc.x(0)
    
    qc = qc.to_gate(label="C")
    return qc 

def V(n):

    qc = QuantumCircuit(n + 1)

    for i in range(n):
        qc.x(i)
    
    qc.mcx([i for i in range(n)], n)

    for i in range(n +1):
        qc.x(i)

    qc = qc.to_gate(label="V")
    return qc



def fsmAlgorithm(x, y, d):

    numIter = floor(log2(d)) + 1
    numDbit = len(bin(d)[2:])
    n = len(x)

    if(n != len(y) or d == 0):
        print("error on input parameter")
        return
    
    qc = initializeCircuit(x, y, d)

    qc.barrier()

    rev = REV(numDbit)
    qc = qc.compose(rev, [i for i in range(2*n, 2*n + numDbit)])
   
    qc.barrier()

    m = M(n)
    qc = qc.compose(m, [i for i in range(0, 2*n)] + [i + 2*n + numDbit for i in range(n)])

    qc. barrier()

    start = 0
    pos = 2*n + numDbit

    for i in range(numIter - 1):
        ext = EXT(n, i)
        qc = qc.compose(ext,[j +  pos + start for j in range(n)] +
                         [j + pos + start + n for j in range(n - 2**i)])
        start+= n 
        n = n - 2**i
        qc.barrier()

    start+=n

    n = len(x)
    initializeD = ID(n)
    qc = qc.compose(initializeD, [pos + start + j for j in range(n)])


    nCurr = n
    currStart = 0

    for i in range(1, numDbit +1):
        
        bcj = BCJ(nCurr, n)
        qc = qc.compose(bcj, [pos - numDbit + i - 1] + [j + pos + currStart for j in range(nCurr)] 
                        + [j + start + pos + (i-1)*n for j in range(n)] + [j + start + pos + i*n for 
                                                                           j in range(n)]) 
        
        rot = ROT(n, 2**(i-1))
        qc = qc.compose(rot, [pos - numDbit + i - 1] + [pos + start + i*n + j 
                                                        for j in range(n)])
        
        c = C(n)
        qc = qc.compose(c, [pos - numDbit + i - 1] +  [j + start + pos + (i-1)*n for j in range(n)] + [j + start + pos + i*n for 
                                                                           j in range(n)])

        currStart+=nCurr
        nCurr = nCurr - 2**i

        '''
        if(i == 2):
            qc.measure([pos + start + i*n + j for j in range(n)], [j for j in range(n)])
            res = u.run(qc, 1)
            u.counts(res)
            return 
        '''

        qc.barrier()
    
    v = V(n)
    qc = qc.compose(v, [j + pos + start + n*(numDbit) for j in range(n+1)])

    qc.measure((pos + start + n*(numDbit + 1)), 0)

    res = u.run(qc, 3)
    u.counts(res)


