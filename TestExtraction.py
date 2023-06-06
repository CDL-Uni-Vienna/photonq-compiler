import pyzx as zx
from extraction.focus_gflow_graph import *
from extraction.perceval import *
import perceval as pcvl
from qiskit import transpile, QuantumCircuit
from qiskit.circuit.library.basis_change import qft

if __name__ == '__main__':
    # circuit='./circuits/grover-orig.qasm'
    circuit = './circuits/grover.qasm'
    # circuit = './circuits/vbe_adder_3.qasm'
    # circuit='QASMBench/small/bell_n4/bell_n4.qasm'
    # circuit='QASMBench/small/linearsolver_n3/linearsolver_n3.qasm'
    ft = qft.QFT(4)
    ft = transpile(ft, basis_gates=['cx', 'cz', 'rx', 'rz', 'h'])

    qc = QuantumCircuit.from_qasm_file(circuit)
    qc = transpile(qc, basis_gates=['cx', 'cz', 'rx', 'rz', 'h'])
    c = zx.Circuit.from_qasm(qc.qasm())
    #c = zx.Circuit.from_qasm(ft.qasm())
    g = c.to_graph().copy()
    zx.draw(g, labels=True).show()

    p = PercevalExtraction(g.to_json())
    zx.draw(p.graph, labels=True).show()
    print(f'Graph: inputs: {p.graph.inputs()}, outputs {p.graph.outputs()}')
    p.extract_clusters_from_graph_ghz_first()
    for c in p.experiment.get_all_clusters().values():
        if c.graph is not None:
            zx.draw(c.graph, labels=True).show()
    exp = p.create_setup(merge=True)
    print("GHZ-Clusters")
    print(p.experiment.ghz_circuits)
    print("Lin Clusters")
    print(p.experiment.lin_circuits)
    print("Fusion: " + str(len(p.experiment.fusion_circuits)))
    print(p.experiment.fusion_circuits)
    # print(p.fusion_circuits)
    ca, exp = p.run()
    print("Lets run the job....")
    pcvl.pdisplay(exp)
    #ca.compute()
    #pcvl.pdisplay(ca.distribution())
    print("Done")
    # print("Photons: " + str(exp.m))
    # print("Optical Components: " + str(exp.ncomponents()))