import sys
import pyzx as zx
from pyzx.rules import remove_ids, spider, lcomp, pivot, MatchLcompType, MatchSpiderType, MatchIdType, MatchPivotType
from pyzx.simplify import simp, Stats, to_gh
from pyzx.graph.base import BaseGraph, VT, ET
from pyzx.utils import VertexType, EdgeType
from pyzx.gflow import gflow
from typing import Callable, Optional, List
from fractions import Fraction

def z_vertices(g: BaseGraph[VT,ET]):
    return g.vertex_set().difference(g.inputs()).difference(g.outputs())

def cluster_clifford_simp(g: BaseGraph[VT,ET], quiet:bool=False, stats:Optional[Stats]=None) -> int:
    spider_simp_4ary(g, quiet=quiet, stats=stats)
    spider_simp_4ary_ext(g, quiet=quiet, stats=stats) #TODO: No proof that this preserves gflow
    to_gh(g)
    i = 0
    while True:
        i1 = id_simp_4ary(g, quiet=quiet, stats=stats)
        i2 = spider_simp_4ary(g, quiet=quiet, stats=stats)
        spider_simp_4ary_ext(g, quiet=quiet, stats=stats) #TODO: Check if necessary (id_simp_4ary may prevent such cases)
        i3 = pivot_simp_4ary(g, quiet=quiet, stats=stats)
        i4 = lcomp_simp_4ary(g, quiet=quiet, stats=stats)
        if i1+i2+i3+i4==0: break
        i += 1
    return i

def lcomp_simp_4ary(g: BaseGraph[VT,ET], matchf:Optional[Callable[[VT],bool]]=None, quiet:bool=False, stats:Optional[Stats]=None) -> int:
    return simp(g, 'lcomp_simp', match_lcomp_4ary_spiders, lcomp, matchf=matchf, quiet=quiet, stats=stats)

def match_lcomp_4ary_spiders(g: BaseGraph[VT,ET]) -> List[MatchLcompType[VT]]:
    candidates = z_vertices(g)
    types = g.types()
    m = []
    while len(candidates) > 0:
        v = candidates.pop()
        va = g.phase(v)

        if not (va == Fraction(1,2) or va == Fraction(3,2)): continue

        if not (
            all(g.edge_type(e) == EdgeType.HADAMARD for e in g.incident_edges(v))
            ): continue
                
        vn = list(g.neighbors(v))

        if not all(types[n] == VertexType.Z for n in vn): continue

        if not check_4ary_lcomp(g, v, vn): continue

        for n in vn: candidates.discard(n)
        m.append((v,vn))
    return m

def check_4ary_lcomp(g: BaseGraph[VT,ET], v: VT, vn: List[VT]) -> bool:
    neighbor_set = set(vn)
    for n in vn:
        unconnected_lc_neighbors = neighbor_set.difference(set(g.neighbors(n))).difference(set([n]))
        non_lc_neighbors = set(g.neighbors(n)).difference(neighbor_set).difference(set([v]))
        if len(unconnected_lc_neighbors) + len(non_lc_neighbors) > 4:
            return False
    return True

def spider_simp_4ary(g: BaseGraph[VT,ET], matchf:Optional[Callable[[VT],bool]]=None, quiet:bool=False, stats:Optional[Stats]=None) -> int:
    return simp(g, 'spider_simp', match_spider_4ary_fusion, spider, matchf=matchf, quiet=quiet, stats=stats)

def match_spider_4ary_fusion(g: BaseGraph[VT,ET]) -> List[MatchSpiderType[VT]]:
    candidates = g.edge_set()
    types = g.types()

    m = []
    while len(candidates) > 0:
        e = candidates.pop()
        if g.edge_type(e) != EdgeType.SIMPLE: continue
        v0, v1 = g.edge_st(e)
        if len(g.neighbors(v0)) + len(g.neighbors(v1)) - 2 > 4: continue
        v0t = types[v0]
        v1t = types[v1]
        if (v0t == v1t and v0t == VertexType.Z):
                for v in g.neighbors(v0):
                    for c in g.incident_edges(v): candidates.discard(c)
                for v in g.neighbors(v1):
                    for c in g.incident_edges(v): candidates.discard(c)
                m.append((v0,v1))
    return m

def id_simp_4ary(g: BaseGraph[VT,ET], matchf:Optional[Callable[[VT],bool]]=None, quiet:bool=False, stats:Optional[Stats]=None) -> int:
    return simp(g, 'id_simp', match_ids_4ary_fusion, remove_ids, matchf=matchf, quiet=quiet, stats=stats)

def match_ids_4ary_fusion(g: BaseGraph[VT,ET]) -> List[MatchIdType[VT]]:

    candidates = g.vertex_set()
    types = g.types()
    phases = g.phases()

    i = 0
    m:List[MatchIdType[VT]] = []

    while len(candidates) > 0:
        v = candidates.pop()
        if phases[v] != 0 or types[v] != VertexType.Z:
            continue
        vn = g.neighbors(v)
        if len(vn) != 2: continue
        v0, v1 = vn
        if len(g.neighbors(v0)) + len(g.neighbors(v1)) - 2 > 4: continue
        if (types[v1] == VertexType.BOUNDARY or types[v0] == VertexType.BOUNDARY):
            continue
        candidates.discard(v0)
        candidates.discard(v1)
        if g.edge_type(g.edge(v,v0)) != g.edge_type(g.edge(v,v1)): #exactly one of them is a hadamard edge
            m.append((v,v0,v1,EdgeType.HADAMARD))
        else: m.append((v,v0,v1,EdgeType.SIMPLE))
        i += 1
    return m

def match_spider_4ary_fusion_ext(g: BaseGraph[VT,ET]):
    res = []
    candidates = filter(lambda e: g.edge_type(e) == EdgeType.SIMPLE, g.edges())
    for candidate in candidates:
        v0, v1 = g.edge_st(candidate)
        if g.types()[v0] == VertexType.Z and g.types()[v1] == VertexType.Z:
            res.append(candidate)
    return res

def insert_spider(g: BaseGraph[VT,ET], matches: List[ET]):
    rem_edges = []
    for m in matches:
        v0, v1 = g.edge_st(m)
        v2 = g.add_vertex(VertexType.Z,g.qubits()[v0],g.rows()[v0]+0.5)
        rem_edges.append(m)
        g.add_edge(g.edge(v0,v2), EdgeType.HADAMARD)
        g.add_edge(g.edge(v2,v1), EdgeType.HADAMARD)
        # etab.append()
    return [{}, [], rem_edges, False]

def spider_simp_4ary_ext(g: BaseGraph[VT,ET], matchf:Optional[Callable[[VT],bool]]=None, quiet:bool=False, stats:Optional[Stats]=None) -> int:
    return simp(g, 'spider_simp', match_spider_4ary_fusion_ext, insert_spider, matchf=matchf, quiet=quiet, stats=stats)

def pivot_simp_4ary(g: BaseGraph[VT,ET], matchf:Optional[Callable[[ET],bool]]=None, quiet:bool=False, stats:Optional[Stats]=None) -> int:
    return simp(g, 'pivot_simp', match_pivot_4ary, pivot, matchf=matchf, quiet=quiet, stats=stats)

def match_pivot_4ary(g: BaseGraph[VT,ET]) -> List[MatchPivotType[VT]]:

    candidates = g.edge_set()
    types = g.types()
    phases = g.phases()
    
    i = 0
    m = []
    while len(candidates) > 0:
        e = candidates.pop()
        if g.edge_type(e) != EdgeType.HADAMARD: continue
        v0, v1 = g.edge_st(e)

        if not (types[v0] == VertexType.Z and types[v1] == VertexType.Z): continue

        v0a = phases[v0]
        v1a = phases[v1]
        if not ((v0a in (0,1)) and (v1a in (0,1))): continue

        if not check_4ary_pivot(g, e): continue

        invalid_edge = False

        v0n = list(g.neighbors(v0))
        v0b = []
        for n in v0n:
            et = g.edge_type(g.edge(v0,n))
            if types[n] == VertexType.Z and et == EdgeType.HADAMARD: pass
            elif types[n] == VertexType.BOUNDARY: v0b.append(n)
            else:
                invalid_edge = True
                break

        if invalid_edge: continue

        v1n = list(g.neighbors(v1))
        v1b = []
        for n in v1n:
            et = g.edge_type(g.edge(v1,n))
            if types[n] == VertexType.Z and et == EdgeType.HADAMARD: pass
            elif types[n] == VertexType.BOUNDARY: v1b.append(n)
            else:
                invalid_edge = True
                break

        if invalid_edge: continue
        if len(v0b) + len(v1b) > 1: continue

        i += 1
        for v in v0n:
            for c in g.incident_edges(v): candidates.discard(c)
        for v in v1n:
            for c in g.incident_edges(v): candidates.discard(c)
        b0 = list(v0b)
        b1 = list(v1b)
        m.append((v0,v1,b0,b1))
    return m

def check_4ary_pivot(g: BaseGraph[VT,ET], e: ET) -> bool:
    v0, v1 = g.edge_st(e)
    n0_set = set(g.neighbors(v0))
    n1_set = set(g.neighbors(v1))
    a = n0_set.difference(n1_set).difference(set([v1]))
    b = n1_set.difference(n0_set).difference(set([v0]))
    c = set(n0_set).intersection(n1_set)
    for s,o in [(a,[b,c,set([v0])]),(b,[a,c,set([v1])]),(c,[a,b,set([v0,v1])])]:
        for v in s:
            unconnected_p_neighbors = o[0].union(o[1]).difference(set(g.neighbors(v))).difference(set([v]))
            non_p_neighbors = set(g.neighbors(v)).difference(o[0].union(o[1])).difference(o[2])
            if len(unconnected_p_neighbors) + len(non_p_neighbors) > 4:
                return False
    return True

def insert_identity(g, v1, v2) -> int:
    orig_type = g.edge_type(g.edge(v1, v2))
    if g.connected(v1, v2):
        g.remove_edge(g.edge(v1, v2))
    vmid = g.add_vertex(VertexType.Z,g.qubits()[v1],g.rows()[v1]+0.5)
    g.add_edge((v1,vmid), EdgeType.HADAMARD)
    if orig_type == EdgeType.HADAMARD:
        g.add_edge((vmid,v2), EdgeType.SIMPLE)
    else:
        g.add_edge((vmid,v2), EdgeType.HADAMARD)
    return vmid

def clean_boundaries(g: BaseGraph[VT,ET]):
    for boundary in g.inputs():
        z_neighbor = list(g.neighbors(boundary))[0]
        e = g.edge(boundary,z_neighbor)
        if g.edge_type(e) == EdgeType.HADAMARD:
            nv = g.add_vertex(VertexType.Z,g.qubits()[boundary],g.rows()[boundary]+0.5)
            g.add_edge(g.edge(boundary, nv), EdgeType.SIMPLE)
            g.add_edge(g.edge(nv, z_neighbor), EdgeType.HADAMARD)
            g.remove_edge(g.edge(boundary, z_neighbor))
    for boundary in g.outputs():
        z_neighbor = list(g.neighbors(boundary))[0]
        e = g.edge(boundary,z_neighbor)
        if g.edge_type(e) == EdgeType.HADAMARD:
            nv = g.add_vertex(VertexType.Z,g.qubits()[boundary],g.rows()[boundary]-0.5)
            g.add_edge(g.edge(boundary, nv), EdgeType.SIMPLE)
            g.add_edge(g.edge(nv, z_neighbor), EdgeType.HADAMARD)
            g.remove_edge(g.edge(boundary, z_neighbor))
    

def to_2d_cluster(g: BaseGraph[VT,ET]):
    gf = gflow(g)
    for v in z_vertices(g):
        vn = g.neighbors(v)
        if len(vn) > 4: #split neighbors according to gflow ordering
            gf_tuples = list(map(lambda n: (n, gf[0][n] if n in gf[0] else -1), vn))
            gf_tuples.sort(key=lambda x: x[1])
            vn1, vn2 = gf_tuples[len(gf_tuples)//2:], gf_tuples[:len(gf_tuples)//2]
            # vn1 = list(filter(lambda n: n in gf[0] and gf[0][n] <= gf[0][v], vn))
            # vn2 = list(filter(lambda n: not n in gf[0] or gf[0][n] > gf[0][v], vn))
            print("split spider: ",v,"in two neighbor sets: ",vn1," and ", vn2)
            split_spider = g.add_vertex(VertexType.Z,g.qubits()[v],g.rows()[v]+0.5,0)
            insert_identity(g, v, split_spider)
            for n in vn2:
                g.remove_edge(g.edge(v, n[0]))
                g.add_edge(g.edge(split_spider, n[0]), EdgeType.HADAMARD)
            gf = gflow(g)
            if not gf:
                print("error: broken gflow")
                import pdb
                pdb.set_trace()

if __name__ == "__main__":
    c = zx.Circuit.load(sys.argv[1]).to_basic_gates().split_phase_gates()
    g = c.to_graph()
    g = zx.simplify.teleport_reduce(g)
    g.track_phases = False
    cluster_clifford_simp(g, quiet=True)
    clean_boundaries(g)
    print(g.to_json())