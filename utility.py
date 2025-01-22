from qiskit import transpile
from qiskit_aer import Aer
import matplotlib.image as mpimg
import matplotlib.pyplot as plt  

def run(qc, shots):
    simulator = Aer.get_backend('aer_simulator')
    compiledCircuit = transpile(qc, simulator)
    job = simulator.run(compiledCircuit, shots = shots).result()
    result = job.get_counts(compiledCircuit)

    print("the circuit depth is : " + str(qc.depth()))
    print("the circuit depth is : " + str(compiledCircuit.depth()))
    return result

def counts(results):
    print("COUNT:")
    for i in results :
        print(" -- " + str(i)[::-1] + " : " + str(results[i]))
    
def draw_circuit(circuit):
    circuit.draw(output='mpl',filename="image.jpg")
    fig, ax = plt.subplots()
    im = mpimg.imread('image.jpg')
    ax.axis('off')
    ax.imshow(im)