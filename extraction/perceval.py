import itertools

import numpy as np
import pyzx as zx
from Experiment import Experiment,PhotonType,CircuitType
import perceval.components as symb
import perceval as pcvl

# makes a ZX-graph "graph-like", such that it only consists of green spiders and hadamards
def to_graph_like(g):
    zx.spider_simp(g)
    zx.to_gh(g)
    zx.id_simp(g)
    zx.spider_simp(g)

# This class is to extract a Perceval Circuit out of a pyzx-graph.
# The graph is described by a json-string, as defined by the pyzx-library.
class PercevalExtraction:

    def __init__(self, graph_json):
        self.total_photons=0
        self.experiment = Experiment()
        self.pcvl_exp = None

        g = zx.Graph.from_json(graph_json)
        g = g.copy()
        to_graph_like(g)
        g = self.optimize_graph(g)
        g.normalize()
        self.graph=g.copy()

    def optimize_graph(self, g):
        zx.simplify.clifford_simp(g)
        zx.full_reduce(g)
        zx.simplify.interior_clifford_simp(g)
        return g

    # Extracts Clusters from the ZX Graph
    # In this function, we first extract all possible GHZ Clusters and Linear Clusters afterwards
    # @max_ghz and @max_lin define the maximum size of ghz and lin clusters respectively
    def extract_clusters_from_graph_ghz_first(self, max_ghz=np.inf, max_lin=np.inf):
        w_graph = self.graph.copy()
        w_graph.remove_vertices(w_graph.inputs() + w_graph.outputs())
        neighbors_list = []
        experiment = Experiment()

        # prepare a list sorted by number of neighbors to get the largest GHZ-Clusters first
        for v in w_graph.vertices():
            if v not in w_graph.inputs() and v not in w_graph.outputs() and len(w_graph.neighbors(v)) > 2:
                neighbors_list.append((v, len(w_graph.neighbors(v))))

        # Find GHZ-Clusters in the Graph und remove their edges from the Graph
        for (v, l) in sorted(neighbors_list, key=lambda tup: tup[1], reverse=True):
            if v not in w_graph.inputs() and v not in w_graph.outputs() and len(w_graph.neighbors(v)) > 2:
                cluster = [v]
                outs = [v] if self.graph.neighbors(v) in self.graph.outputs() else []
                i = 0
                while len(w_graph.neighbors(v)) > 0 and i < max_ghz:
                    n = list(w_graph.neighbors(v))[0]
                    if n not in w_graph.outputs() and n not in w_graph.inputs():
                        cluster.append(n)
                        i = i + 1
                    w_graph.remove_edge(w_graph.edge(v, n))
                    if any(item in self.graph.neighbors(n) for item in self.graph.outputs()):
                        outs.append(n)
                tmp_graph = self.extract_subgraph(cluster)
                experiment.add_cluster(CircuitType.GHZ, graph=tmp_graph, nodes=cluster, readouts=outs)

        # Find Lin-Clusters in the remaining Graph und remove their edges from the Graph.
        while w_graph.num_edges() != 0:
            for v in list(w_graph.vertices()):
                if len(w_graph.neighbors(v)) == 1:
                    n = v
                    cluster = [n]
                    outs = [n] if self.graph.neighbors(n) in self.graph.outputs() else []
                    i = 0
                    while len(w_graph.neighbors(n)) > 0 and i < max_lin:
                        n2 = list(w_graph.neighbors(n))[0]
                        cluster.append(n2)
                        w_graph.remove_edge(w_graph.edge(n, n2))
                        if any(item in self.graph.neighbors(n2) for item in self.graph.outputs()):
                            outs.append(n2)
                        n = n2
                        i = i + 1
                    tmp_graph = self.extract_subgraph(cluster)

                    # we consider everything smaller than 4, a GHZ
                    if len(cluster) == 3:
                        # in case of 3, we have to restructer the cluster to have the central node on the first position
                        experiment.add_cluster(CircuitType.GHZ, graph=tmp_graph,
                                               nodes=[cluster[1],cluster[0], cluster[2]], readouts=outs)
                    elif len(cluster) == 2:  # we consider everything smaller than 4, a GHZ
                        experiment.add_cluster(CircuitType.GHZ, graph=tmp_graph,nodes=cluster, readouts=outs)
                    else:
                        experiment.add_cluster(CircuitType.LIN, graph=tmp_graph, nodes=cluster, readouts=outs)

            # There might be some circles left. if so, extract one single 3-qubit ghz state before continue with linear
            # clusters.
            for v in list(w_graph.vertices()):
                if len(w_graph.neighbors(v)) == 2:
                    neighbors = list(w_graph.neighbors(v))
                    cluster = [v] + neighbors
                    outs = []
                    for n in cluster:
                        if any(item in self.graph.neighbors(n)  for item in self.graph.outputs()):
                            outs.append(n)
                    w_graph.remove_edge(w_graph.edge(v, neighbors[0]))
                    w_graph.remove_edge(w_graph.edge(v, neighbors[1]))

                    tmp_graph = self.extract_subgraph(cluster)
                    experiment.add_cluster(CircuitType.GHZ, graph=tmp_graph, nodes=cluster, readouts=outs)
                    break
        experiment.add_fusions()
        self.experiment = experiment

    # Creates a Perceval Setup for the Experiment
    def create_setup(self, merge=False):

        self.total_photons = len(self.experiment.photons.items())
        exp = symb.Circuit(self.total_photons, name="CDL Experiment")

        # Create circuits for all the GHZ clusters
        for circ in self.experiment.ghz_circuits.values():
            exp.add(circ.photons[0].exp_id, circ.circuit, merge)


        # Create circuits for all the Lin clusters
        for circ in self.experiment.lin_circuits.values():
            exp.add(circ.photons[0].exp_id, circ.circuit, merge)

        exp.add(0, symb.PERM(list(range(self.total_photons))))
        # Create circuits for all the Fusions
        for circ in self.experiment.fusion_circuits.values():
            exp.add(circ.photons[0].exp_id, circ.circuit, merge)

        exp.add(0,symb.PERM(list(range(self.total_photons)))) ## just a barrier
        # add the projective measurements
        for phot in self.experiment.photons.values():
            if phot.type in (PhotonType.COMP, PhotonType.READOUT):
                exp.add(phot.exp_id, symb.QWP(phot.angle[0]))
                exp.add(phot.exp_id, symb.HWP(phot.angle[1]))

        self.pcvl_exp = exp

        return exp

    def extract_subgraph(self, cluster_list, head=-1):
        # if it has a head, it's a ghz state
        if head != -1:
            edges = [(head, i) for i in cluster_list if i != head] + [(i, head) for i in cluster_list if i != head]
        else:
            edges = [(cluster_list[i],cluster_list[i+1]) for i in range(len(cluster_list)) if i+1 < len(cluster_list)] + \
                    [(cluster_list[i+1],cluster_list[i]) for i in range(len(cluster_list)) if i+1 < len(cluster_list)]
        tmp_graph = self.graph.copy()
        tmp_graph.remove_edges([i for i in list(tmp_graph.edges()) if i not in edges])
        tmp_graph.remove_vertices([i for i in list(tmp_graph.vertices()) if i not in cluster_list])
        return tmp_graph


    # create the input state. for each photon of the original circuit, the input state is defined by the "Photon"
    # object (typically |{P:H}>. For the ancillary ones for measurement the input is "|0>"
    def input_states(self):
        in_state = "|"
        in_state = in_state + ",".join([f'{ph.in_state}' for ph in self.experiment.photons.values()])
        in_state = in_state + "," + ",".join([f'{item}' for item in [0] * self.total_photons])
        in_state = in_state + ">"

        return [pcvl.BasicState(in_state)]

    # Create the output states.
    # For each WITNESS photon, the output is fixed to "|{P:H}>" and translated
    # to |0,0> (as we have to measure two photons).
    # For Computings Photons, the output is fixed to |P:V>, aka |1,0>
    # For the READOUTS, the output can be either "|{P:H}>" or "|{P:V}>
    def output_states(self):
        (readouts, witnesses, comps) = (self.experiment.get_readouts(), self.experiment.get_witnesses(),
                                        self.experiment.get_computes())
        out_states = {}
        x = 0
        (zero, one) = ([0, 1], [1, 0])
        for st in itertools.product([zero, one], repeat=len(readouts)):
            basic_out_state = [[]] * self.total_photons
            for w in witnesses:
                basic_out_state[w.exp_id] = zero
            for c in comps:
                basic_out_state[c.exp_id] = zero
            for i in range(len(readouts)):
                basic_out_state[readouts[i].exp_id] = st[i]
            out_states[
                pcvl.BasicState(list(itertools.chain.from_iterable(basic_out_state)))] = f'|{x:0{len(readouts)}b}>'
            x = x + 1
        return out_states

    def run(self):
        if self.pcvl_exp is None or self.experiment is None:
            print("No experiment ot run")
            return

        # We have to increase the size by 2, as Perceval currently does not support measurement in polarization
        # hence, we add a second rail for each photon and a PBS before measurement
        # we add the rails at the end of the circuit and permutate them upwards just before the measurement,
        # e.g. in case of 10 photons, photon 10 will be the helper for photon 0; photon 11 for 1 and so on
        phs = self.total_photons
        exp = symb.Circuit(phs*2, name="CDL Experiment")
        exp.add(0,self.pcvl_exp)
        exp.add(0,symb.PERM(sum(list([2*i] for i in range(0,phs)) + list([2*i+1] for i in range(0,phs)), [])))
        for i in range(0,phs):
            exp.add(i*2, symb.PBS())

        sim = pcvl.Processor("Naive", exp)
        ca = pcvl.algorithm.Analyzer(sim,input_states=self.input_states(),output_states=self.output_states())
        return ca, exp



