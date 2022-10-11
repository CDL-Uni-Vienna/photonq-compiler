import os, sys
import pyzx as zx
import perceval.lib.phys as phys
from qiskit import QuantumCircuit, transpile
from extraction.perceval import PercevalExtraction
from perceval.converters import QiskitConverter
from optimization.optimizer import cluster_clifford_simp, clean_boundaries, spider_simp_4ary

small_circuits = ['wstate_n3','adder_n4','cat_state_n4','deutsch_n2','fredkin_n3','grover_n2','hs4_n4','inverseqft_n4','iswap_n2','linearsolver_n3','qft_n4','teleportation_n3','toffoli_n3','shor_n5','bell_n4']
medium_circuits = ['adder_n10','bv_n14','cc_n12','multiply_n13','qf21_n15','qft_n15','qpe_n9','sat_n11','simon_n6','multiplier_n15']

def evaluate_folder(source_folder,target_folder):
    list_of_files = filter( lambda x: os.path.isdir(os.path.join(source_folder, x)),os.listdir(source_folder) )
    for f in list_of_files:
        if not f in small_circuits and not f in medium_circuits:
            # print("skip ",f)
            continue
        full_path = os.path.join(source_folder,f,f+".qasm")
        # print("evaluate ", full_path)
        try:
            qc = QuantumCircuit.from_qasm_file(full_path)
        except:
            print("qiskit failed to parse this circuit")
            continue
        try:
            qct = transpile(qc, basis_gates=['cx', 'cz', 'rz', 'rx', 'h'])
        except:
            print("qiskit failed to transpile this circuit")
            continue            
        qasm_str = qct.qasm()
        try:
            c = zx.Circuit.from_qasm(qasm_str)
        except:
            # print("pyzx failed to parse this circuit")
            continue
        g = c.to_graph()
        # spider_simp_4ary(g, quiet=True)
        # zx.simplify.to_gh(g)
        # print("before ", g)
        g = zx.simplify.teleport_reduce(g)
        g.track_phases = False
        zx.simplify.clifford_simp(g)
        # cluster_clifford_simp(g, quiet=True)
        # clean_boundaries(g)
        p = PercevalExtraction(g.to_json())
        p.extract_clusters_from_graph_ghz_first()
        
        # print("GHZ-Clusters: " + str(len(p.ghz_clusters_dict)))
        # print(p.ghz_clusters_dict)
        # print("Lin Clusters: " + str(len(p.lin_clusters_dict)))
        # print(p.qubit_fusions)
        # print("Fusion: " + str(len(p.qubit_fusions)))
        # print(p.qubit_fusions)
        exp, clusters = p.create_setup(True, False)
        # print("Photons: " + str(exp.m))
        # print("Optical Components: " + str(exp.ncomponents()))
        # with open(target_folder+f+"out.json", "w") as outfile:
        #     outfile.write(g.to_json())
        # print("after ", g)

        ## Perceval Extraction ####
        qct = transpile(qc, basis_gates=['cx', 'rx', 'rz', 'h'])
        qct.remove_final_measurements()
        qiskit_converter = QiskitConverter(phys)
        try:
            qp = qiskit_converter.convert(qct, heralded=False)
            qp_comps, qp_phot = qp.circuit.ncomponents(), qp.circuit.m
        except:
            qp_comps, qp_phot = "--", "--"
        try:
            herald = qiskit_converter.convert(qct, heralded=True)
            her_comps, her_phot = herald.circuit.ncomponents(), herald.circuit.m
        except:
            her_comps, her_phot = "--", "--"


        ## Result as Latex Table

        print(f.replace("_",'-')," & ",len(p.ghz_clusters_dict)," & ",len(p.lin_clusters_dict)," & ",exp.m," & ",
              exp.ncomponents(), " & ", qp_phot, " & ", qp_comps," & ", her_phot, " & ", her_comps , "\\\\")



if __name__ == "__main__":
    #evaluate_folder(sys.argv[1], sys.argv[2])
    evaluate_folder('D:/Tobias/Dokumente/workspaces/photonq-compiler/QASMBench/medium/', '')