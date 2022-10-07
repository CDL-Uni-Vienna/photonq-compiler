import numpy as np
import pyzx as zx
import perceval as pcvl
import perceval.lib
import perceval.lib.phys as symb
import sympy as sp
from .focus_gflow_graph import *

# Returns a Perceval Circuit for fusing two photons.
# If the photons are not neighbors, we have to swap (permutate) the ph2 with the photon next to ph1, do
# the fusion and swap back.
def fusion(ph1,ph2):
    l = abs(ph2-ph1)
    circ = pcvl.Circuit(m=l+1, name="Fuse " + str(ph1) + "-" + str(ph2))
    if l > 1:
        a,*b,c = list(range(0,l))
        perm=symb.PERM([c,*b,a])
        circ.add(1,perm)
    circ.add((0,1), symb.PBS())
    circ.add(1, symb.HWP(sp.pi / 8)) #TODO remove if we wanna measure in A/D basis
    if l > 1:
        circ.add(1,perm)
    return circ

# Retuns a Perceval Circuit for a Linear cluster with a "name" and the given number of photons
def linear_cluster(name, photons):
    circ = symb.Circuit(m=photons, name="LIN " + str(name))
    for i in range(photons):
        circ.add(i, symb.HWP(sp.pi / 8))

    for i in range(photons-1):
        circ.add((i, i + 1), symb.PBS())
        circ.add(i+1, symb.HWP(sp.pi / 8))
    circ.add(0, symb.HWP(sp.pi / 8))


    return circ

# Retuns a Perceval Circuit for a GHZ cluster with a "name" and the given number of photons
def ghz_cluster(name, photons):
    circ = symb.Circuit(m=photons, name="GHZ " + str(name))
    for i in range(photons):
        circ.add(i, symb.HWP(sp.pi / 4))

    for i in range(photons-1):
        circ.add((i, i + 1), symb.PBS())

    for i in range(photons-1):
        circ.add(i+1, symb.HWP(sp.pi / 8))

    return circ

# makes a ZX-graph "graph-like", such that it only consists of green spiders and hadamards
def to_graph_like(g):
    zx.spider_simp(g)
    zx.to_gh(g)
    zx.id_simp(g)
    zx.spider_simp(g)

# This class is to extract a Perceval Circuit out of a pyzx-graph.
# The graph is described by a json-string, as defined by the pyzx-library.
class PercevalExtraction:
    lin_clusters=[]
    ghz_clusters=[]
    lin_clusters_dict={}
    ghz_clusters_dict={}
    horizontal_fusions=[]
    qubit_fusions = []
    total_photons=0
    photonic_clusters=[]
    graph = zx.Graph()
    ghz_circuits = {}
    lin_circuits = {}
    fusion_circuits = {}
    experiment = None

    def __init__(self, graph_json):
        g = zx.Graph.from_json(graph_json)
        to_graph_like(g)
        g = g.copy()
        zx.simplify.interior_clifford_simp(g)
        self.graph=g

    # Extracts Clusters from the ZX Graph
    # In this function, we first extract all possible GHZ Clusters and Linear Clusters afterwards
    # @max_ghz and @max_lin define the maximum size of ghz and lin clusters respectively
    def extract_clusters_from_graph_ghz_first(self, max_ghz = np.inf, max_lin = np.inf):
        w_graph = self.graph.copy()
        w_graph.remove_vertices(w_graph.inputs() + w_graph.outputs())
        neighbors_list = []
        lin_clusters = {}
        ghz_clusters = {}
        fusions = {}

        # prepare a list sorted by number of neighbors to get the largest GHZ-Clusters first
        for v in w_graph.vertices():
            if v not in w_graph.inputs() and v not in w_graph.outputs() and len(w_graph.neighbors(v)) > 2:
                neighbors_list.append((v, len(w_graph.neighbors(v))))

        # Find GHZ-Clusters in the Graph und remove their edges from the Graph
        for (v, l) in sorted(neighbors_list, key=lambda tup: tup[1], reverse=True):
            if v not in w_graph.inputs() and v not in w_graph.outputs() and len(w_graph.neighbors(v)) > 2:
                ghz_clusters[v] = [v]
                i = 0
                while len(w_graph.neighbors(v)) > 0 and i < max_ghz:
                    n = list(w_graph.neighbors(v))[0]
                    if n not in w_graph.outputs() and n not in w_graph.inputs():
                        ghz_clusters[v].append(n)
                        i = i + 1
                    w_graph.remove_edge(w_graph.edge(v, n))

        # Find Lin-Clusters in the remaining Graph und remove their edges from the Graph.
        for v in list(w_graph.vertices()):
            if len(w_graph.neighbors(v)) == 1:
                n = v
                lin_clusters[v] = [n]
                i = 0
                while len(w_graph.neighbors(n)) > 0 and i < max_lin:
                    n2 = list(w_graph.neighbors(n))[0]
                    lin_clusters[v].append(n2)
                    w_graph.remove_edge(w_graph.edge(n, n2))
                    n = n2
                    i=i+1

        # There should be no edges left in the Graph....
        if (w_graph.num_edges() == 0):
            print("Everything worked out fine, no. edges = " + str(w_graph.num_edges()))
        else:
            print("Some edges are left..., no. edges = " + str(w_graph.num_edges()))

        # find edges that are connected with two clusters. They have to be "fused"
        # with the extraction method above, we only have to find GHZ->GHZ and GHZ->Lin, because we cannot have two
        # linear cluster connected (if there would be, there is a 'small' GHZ in between them)
        for k, v in ghz_clusters.items():
            # first look for other GHZ clusters
            for k2, v2 in ghz_clusters.items():
                if k == k2:
                    continue
                for i in v:
                    if i in v2:
                        if i not in fusions.keys():
                            fusions[i] = []
                        fusions[i].append(k)
                        fusions[i].append(k2)

            # next look for connected linear clusters
            for k2, v2 in lin_clusters.items():
                if k == k2:
                    continue
                for i in v:
                    if i in v2:
                        if i not in fusions.keys():
                            fusions[i] = []
                        fusions[i].append(k)
                        fusions[i].append(k2)
        for k, v in fusions.items():
            fusions[k] = list(set(v))

        self.qubit_fusions=fusions
        self.ghz_clusters_dict=ghz_clusters
        self.lin_clusters_dict=lin_clusters

    # this method is old and not working at the moment.
    # The idea here is to only create linear clusters and connect them with
    # fusions wherever a vertical connection appears.
    def extract_clusters_from_graph_linears_first(self):
        g2 = self.graph.copy()
        fg_graph = build_focused_gflow_graph(g2, focus_gflow(g2, gflow(g2)))
        i = 0
        lin_clusters = {}
        while len(list(fg_graph.edges())) > 0:
            n = max(fg_graph.vertices())
            lin_clusters[i] = [n]
            while n not in fg_graph.inputs() and len(fg_graph.neighbors(n)) > 0:
                neigh = list(fg_graph.neighbors(n))[0]
                lin_clusters[i].insert(0, neigh)
                fg_graph.remove_edge(fg_graph.edge(n, neigh))
                if len(fg_graph.neighbors(n)) == 0:
                    fg_graph.remove_vertex(n)
                n = neigh
            if len(fg_graph.neighbors(n)) == 0:
                fg_graph.remove_vertex(n)
            i = i + 1

        vertical_fusions = []
        qubit_fusions = []
        for key in lin_clusters.keys():
            for node in lin_clusters[key]:
                for key2 in lin_clusters.keys():
                    if key2 > key:
                        for node2 in lin_clusters[key2]:
                            if g2.connected(node, node2):
                                vertical_fusions.append([(key, node), (key2, node2)])
                        if node in lin_clusters[key2]:
                            qubit_fusions.append([(key, node), (key2, node)])

        self.lin_clusters=lin_clusters
        self.vertical_fusions=vertical_fusions
        self.qubit_fusions=qubit_fusions

    #Creates a Perceval Setup for the graph
    def create_setup(self, merge=False, reorder_before_measurement=True):
        ghz_circuits = {}
        lin_circuits = {}
        fusion_circuits = {}
        clusters={}
        # we need as many photons as items in all GHZ+Lin Clusters
        photons = sum(len(item) for item in self.ghz_clusters_dict.values()) + sum(
            len(item) for item in self.lin_clusters_dict.values())
        self.total_photons = photons
        exp = symb.Circuit(photons, name="CDL Experiment")
        i = 0

        # Create circuits for all the GHZ clusters
        for k, v in self.ghz_clusters_dict.items():
            ghz_circuits[k] = ghz_cluster(k, len(v))
            exp.add(i, ghz_circuits[k], merge)
            tmp = {}
            for j in range(len(v)):
                # "i" is the actual photon number. we need that for the fusion later
                tmp[v[j]] = i
                i = i + 1
            self.ghz_clusters_dict[k]= tmp
            clusters[k] = tmp

        #Create circuits for all the Lin-Clusters
        for k, v in self.lin_clusters_dict.items():
            lin_circuits[k] = linear_cluster(k, len(v))
            exp.add(i, lin_circuits[k], merge)
            tmp = {}
            for j in range(len(v)):
                # "i" is the actual photon number. we need that for the fusion later
                tmp[v[j]] = i
                i = i + 1
            self.lin_clusters_dict[k] = tmp
            clusters[k] = tmp

        #create fusion circuits.
        for k, v in self.qubit_fusions.items():
            fusion_circuits[k] = {"orig": 0, "whitness": [], "circuits": []}
            for i in range(len(v) - 1):
                if i == 0:
                    fusion_circuits[k]["orig"] = clusters[v[i]][k]
                fusion_circuits[k]["whitness"].append(clusters[v[i + 1]][k])
                ph1 = clusters[v[i]][k]
                ph2 = clusters[v[i + 1]][k]
                fuse = fusion(ph1, ph2)
                fusion_circuits[k]["circuits"].append(fuse)
                i = i + 1
                exp.add(min(ph1, ph2), fuse, merge)

        self.ghz_circuits = ghz_circuits
        self.lin_circuits = lin_circuits
        self.fusion_circuits = fusion_circuits

        if reorder_before_measurement:
            perm = self.reorder_before_measurement()
            exp.add(perm, 0)


        self.experiment = exp

        return exp, clusters

    def reorder_before_measurement(self):
        perm = list(range(0, self.total_photons))

        i = self.total_photons-1
        # put all fusions whitnesses to the bottom of the circuit
        for ver in self.fusion_circuits.keys():
            for w in self.get_whitnesses(ver):
                perm[w] = i
                self.update_whitnesses(ver, w, i)
                self.update_photon(w,i)
                perm[i] = w
                self.update_photon(i, w)
                i = i - 1

        # put all outputs aka readouts at the top of the circuit
        i = 0
        for outs in map (lambda o: list(self.graph.neighbors(o))[0], list(self.graph.outputs())):
            if outs in self.fusion_circuits.keys():
                orig = self.update_origin(outs, i)
                perm[i] = orig
                self.update_photon(orig, i)
                perm[orig] = i
                self.update_photon(i, orig)
            else:
                ph = self.get_photon_for_node(outs)
                perm[i] = ph
                self.update_photon(ph, i)
                perm[ph] = i
                self.update_photon(i, ph)
            i = i + 1
        print(perm)
        return symb.PERM(perm)

    def get_whitnesses(self, node_id):
        return self.fusion_circuits[node_id]["whitness"]


    def update_whitnesses(self, node_id, old, new):
        for i in range(len(self.fusion_circuits[node_id]["whitness"])):
            if self.fusion_circuits[node_id]["whitness"][i] == old:
                self.fusion_circuits[node_id]["whitness"][i] = new
                return

    def update_origin(self, node_id, new):
        old = self.fusion_circuits[node_id]["orig"]
        self.fusion_circuits[node_id]["orig"] = new
        return old

    def update_photon(self, old, new):
        for k,v in self.ghz_clusters_dict.items():
            for k2,v2 in v.items():
                if v2 == old:
                    v[k2] = new
                    return

        for k,v in self.lin_clusters_dict.items():
            for k2,v2 in v.items():
                if v2 == old:
                    v[k2] = new

    def get_photon_for_node(self, node_id):
        if node_id in self.fusion_circuits.keys():
            return self.fusion_circuits[node_id]["orig"]
        for k,v in self.ghz_clusters_dict.keys():
            if node_id in v.keys():
                return v[node_id]
        for k, v in self.lin_clusters_dict.keys():
            if node_id in v.keys():
                return v[node_id]

    def run(self):
        if self.experiment is None:
            print("No experiment ot run")

