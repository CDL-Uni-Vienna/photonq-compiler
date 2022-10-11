import pyzx as zx
from extraction.focus_gflow_graph import *
from extraction.perceval import *
import perceval as pcvl
from qiskit import transpile, QuantumCircuit
from qiskit.circuit.library.basis_change import qft

if __name__ == '__main__':
    # circuit='./circuits/grover-orig.qasm'
    circuit = './circuits/grover.qasm'
    circuit = './circuits/vbe_adder_3.qasm'
    circuit='QASMBench/small/bell_n4/bell_n4.qasm'
    ft = qft.QFT(4)
    ft = transpile(ft, basis_gates=['cx', 'cz', 'rx', 'rz', 'h'])

    qc = QuantumCircuit.from_qasm_file(circuit)
    qc = transpile(qc, basis_gates=['cx', 'cz', 'rx', 'rz', 'h'])
    c = zx.Circuit.from_qasm(qc.qasm())
    #c = zx.Circuit.from_qasm(ft.qasm())
    g = c.to_graph()
    #zx.draw(g)

    g1 = g.copy()
    to_graph_like(g1)
    g1 = g1.copy()
    zx.simplify.clifford_simp(g1)
    #zx.full_reduce(g1)
    # g.normalize()
    g1 = g1.copy()
    #zx.draw(zx.Graph.from_json(g1.to_json()), labels=True)

    p = PercevalExtraction(g1.to_json())
    p.extract_clusters_from_graph_ghz_first()
    exp, clusters = p.create_setup(True, False)
    print("GHZ-Clusters: " + str(len(p.ghz_clusters_dict)))
    print(p.ghz_clusters_dict)
    print("Lin Clusters: " + str(len(p.lin_clusters_dict)))
    print(p.lin_clusters_dict)
    print("Fusion: " + str(len(p.qubit_fusions)))
    print(p.qubit_fusions)
    print(p.fusion_circuits)
    # pcvl.pdisplay(exp)
    print("Photons: " + str(exp.m))
    print("Optical Components: " + str(exp.ncomponents()))