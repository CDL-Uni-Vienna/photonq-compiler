import os, sys
import pyzx as zx
from qiskit import QuantumCircuit, transpile
from optimizer import cluster_clifford_simp, clean_boundaries

def evaluate_folder(source_folder,target_folder):
    list_of_files = filter( lambda x: os.path.isdir(os.path.join(source_folder, x)),os.listdir(source_folder) )
    for f in list_of_files:
        full_path = os.path.join(source_folder,f,f+".qasm")
        print("evaluate ", full_path)
        try:
            qc = QuantumCircuit.from_qasm_file(full_path)
        except:
            print("qiskit failed to parse this circuit")
            continue
        try:
            qct = transpile(qc, basis_gates=['cx', 'rz', 'rx', 'h'])
        except:
            print("qiskit failed to transpile this circuit")
            continue            
        qasm_str = qct.qasm()
        try:
            c = zx.Circuit.from_qasm(qasm_str)
        except:
            print("pyzx failed to parse this circuit")
            continue
        g = c.to_graph()
        zx.simplify.spider_simp(g, quiet=True)
        zx.simplify.to_gh(g)
        print("before ", g)
        g = zx.simplify.teleport_reduce(g)
        g.track_phases = False
        cluster_clifford_simp(g, quiet=True)
        clean_boundaries(g)
        with open(target_folder+f+"out.json", "w") as outfile:
            outfile.write(g.to_json())
        print("after ", g)

if __name__ == "__main__":
    evaluate_folder(sys.argv[1],sys.argv[2])