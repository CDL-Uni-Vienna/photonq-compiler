import itertools
from enum import Enum

import pyzx as zx
import perceval as pcvl
import perceval.components as symb
import sympy as sp

class PhotonType(Enum):
    READOUT = "READOUT"
    COMP = "COMP"
    WITNESS = "WITNESS"
    LOSS = "LOSS"
    NONE = "99"

    def __str__(self):
        return self.name


class CircuitType(Enum):
    GHZ = "GHZ"
    LIN = "LIN"
    FUSION = "FUSION"
    NONE = "99"

    def __str__(self):
        return self.name


class Experiment(object):

    def __init__(self):
        self.ghz_circuits = {}
        self.lin_circuits = {}
        self.fusion_circuits = {}
        self.photons = {}
        Photon.newid = itertools.count()
        Cluster.c_newid = itertools.count()

    def add_cluster(self, type:CircuitType, graph:zx.graph.graph.BaseGraph, nodes:list, readouts:list):
        tmp = []
        head = nodes[0]
        for i in nodes:
            if i in readouts:
                ph = Photon(type=PhotonType.READOUT, zx_id=i, angle=graph.phase(vertex=i))
            else:
                ph = Photon(type=PhotonType.COMP, zx_id=i, angle=graph.phase(vertex=i))
            self.photons[ph.id] = ph
            tmp.append(ph)
            if i == head:
                head = ph

        if type ==  CircuitType.GHZ:
            self.ghz_circuits[head] = Cluster(type=type,name=head, graph=graph, photons=tmp)
        elif type == CircuitType.LIN:
            self.lin_circuits[head] = Cluster(type=type,name=head, graph=graph, photons=tmp)
        else:
            print("Type not supported")

    def get_readouts(self):
        return [ph for ph in self.photons.values() if ph.type == PhotonType.READOUT]

    def get_computes(self):
        return [ph for ph in self.photons.values() if ph.type == PhotonType.COMP]

    def get_witnesses(self):
        return [ph for ph in self.photons.values() if ph.type == PhotonType.WITNESS]

    # find edges that are connected with two clusters. They have to be "fused"
    def add_fusions(self):
        for k1, v1 in self.ghz_circuits.items():
            for k2, v2 in self.ghz_circuits.items():
                for ph in Cluster.get_fusion_photons(v1,v2):
                    if (ph[0],ph[1]) not in self.fusion_circuits.keys() and (ph[1],ph[0]) not in self.fusion_circuits.keys():
                        self.fusion_circuits[(ph[0],ph[1])] = Cluster(type=CircuitType.FUSION,name=ph[0].zx_id, photons=ph)

            for k2, v2 in self.lin_circuits.items():
                for ph in Cluster.get_fusion_photons(v1,v2):
                    if (ph[0],ph[1]) not in self.fusion_circuits.keys() and (ph[1],ph[0]) not in self.fusion_circuits.keys():
                        self.fusion_circuits[(ph[0],ph[1])] = Cluster(type=CircuitType.FUSION,name=ph[0].zx_id, photons=ph)

    def get_all_clusters(self):
        return self.ghz_circuits | self.lin_circuits

    def get_all_circuits(self):
        return self.get_all_clusters() | self.fusion_circuits


class Photon(object):
    newid = itertools.count()

    def __init__(self, type:PhotonType, zx_id=0, angle=0):
        self.id = next(Photon.newid)
        self.exp_id = self.id
        self.type = type
        self.zx_id = zx_id
        self.in_state = "{P:H}"
        # angle for QWP and HWP
        if angle:
            self.angle=[sp.pi/4, (sp.pi - (2*angle*sp.pi)) / 8 ]
        else :
            self.angle=[sp.pi/4,sp.pi/8]


    @staticmethod
    def get_photon_ids(photons:list) ->list:
        ret = []
        for ph in photons:
            ret.append(ph.exp_id)
        return ret

    @staticmethod
    def get_node_ids(photons: list) -> list:
        ret = []
        for ph in photons:
            ret.append(ph.zx_id)
        return ret

    @staticmethod
    def get_photon_zx_id(zx_id, photons: list):
        for ph in photons:
            if ph.zx_id == zx_id:
                return ph
        return None

    def __str__(self):
        return f'{str(self.zx_id)}-{str(self.exp_id)}({str(self.type)})'

    def __repr__(self):
        return f'{str(self.zx_id)}-{str(self.exp_id)}({str(self.type)})'

class Cluster(object):
    c_newid = itertools.count()

    def __init__(self, type:CircuitType, name="", graph=None, photons=[]):
        self.id = next(Cluster.c_newid)
        self.name = name
        self.photons = photons
        self.graph = graph
        self.type = type
        if type ==  CircuitType.GHZ:
            self.circuit = ghz_cluster(name, len(photons))
        elif type ==  CircuitType.LIN:
            self.circuit = linear_cluster(name, len(photons))
        elif type ==  CircuitType.FUSION:
            self.circuit = fusion(photons[0], photons[1])
        else:
            self.circuit = None
    def get_photon_zx_id(self, zx_id):
        for ph in self.photons:
            if ph.zx_id == zx_id:
                return ph
        return None

    @staticmethod
    def get_fusion_photons(c1, c2)->list:
        if c1 == c2:
            return []
        ret_list=[]
        for p1 in c1.photons:
            for p2 in c2.photons:
                ph = c2.get_photon_zx_id(p1.zx_id)
                if ph:
                    ret_list.append((p1,ph))
        return ret_list

    def __str__(self):
        return str(self.type) + str(self.photons)


    def __repr__(self):
        return str(self.type) + str(self.photons)

# Returns a Perceval Circuit for fusing two photons.
# If the photons are not neighbors, we have to swap (permutate) the ph2 with the photon next to ph1, do
# the fusion and swap back.
def fusion(ph1:Photon,ph2:Photon):

    l = abs(ph2.exp_id-ph1.exp_id)
    circ = pcvl.Circuit(m=l+1, name="Fuse " + str(ph1.exp_id) + "-" + str(ph2.exp_id))
    if l > 1:
        a,*b,c = list(range(0,l))
        perm=symb.PERM([c,*b,a])
        circ.add(1,perm)
    circ.add((0,1), symb.PBS())
    circ.add(1, symb.HWP(sp.pi / 8)) #TODO remove if we wanna measure in A/D basis
    if l > 1:
        circ.add(1,perm)
    ph2.type=PhotonType.WITNESS
    return circ

# Retuns a Perceval Circuit for a Linear cluster with a "name" and the given number of photons
def linear_cluster(name, photons):
    circ = symb.Circuit(m=photons, name="LIN " + str(name))
    for i in range(photons):
        circ.add(i, symb.HWP(sp.pi / 8))

    for i in range(photons-1):
        circ.add((i, i + 1), symb.PBS())
        if i >=1 and i != photons-2:
            circ.add(i+1, symb.HWP(sp.pi / 8))
    circ.add(0,symb.PERM(list(range(photons)))) # just to have the circuit nicer displayed (like a barrier)
    circ.add(0, symb.HWP(sp.pi / 8))
    circ.add(photons-1,symb.HWP(sp.pi/8))

    return circ

# Retuns a Perceval Circuit for a GHZ cluster with a "name" and the given number of photons
def ghz_cluster(name, photons):
    circ = pcvl.Circuit(m=photons, name="GHZ " + str(name))
    for i in range(photons):
        circ.add(i, symb.HWP(sp.pi / 8))

    for i in range(photons-1):
        circ.add((i, i + 1), symb.PBS())

    # The Hadamard correction before measurement
    for i in range(1,photons):
        circ.add(i, symb.HWP(sp.pi / 8))

    return circ