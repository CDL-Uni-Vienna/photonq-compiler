from typing import List, Dict, Set, Tuple, Optional
from pyzx.graph.base import BaseGraph, VT, ET
from pyzx.graph.graph_s import GraphS
from pyzx.linalg import Mat2, CNOTMaker

'''builds focused gflow graph where vertices are connected iff
 they are connected in the original graph and appear in one of each correction sets'''
def build_focused_gflow_graph(g, gff):
  data = []
  for out in g.outputs():
    line = build_line_rec(g,out, gff)
    data.append(line)
  return build_diagram_from_lines(g, data)

'''focuses a maximally delayed gflow, i.e. the odd neighborhood of the correction set of each vertex is the vertex itself + output vertices'''
def focus_gflow(g: BaseGraph[VT,ET], gf: List[Tuple[Dict[VT,int], Dict[VT,Set[VT]], int]]):
  l:     Dict[VT,int]      = {}
  cset: Dict[VT, Set[VT]] = {}
  for v in g.outputs():
    l[v] = 0
    cset[v] = set() #usually outputs do not occur in g function but we need this construct as a temporary helper
  processed = set(g.outputs())

  unprocessed = set(g.vertices()).difference(processed)
  k = 1
  while True:
    candidates = []
    for v in unprocessed:
      try:
        odd_n = set(get_odd_neighbourhood(g,gf[1][v]))
      except:
        import pdb
        pdb.set_trace()
      odd_n.discard(v)
      no_candidate = False
      for n in odd_n:
        if not n in cset:
          no_candidate = True
          break
      if no_candidate:
        continue

      c_set = set(gf[1][v])
      for n in odd_n:
        c_set.symmetric_difference_update(cset[n])
      cset[v] = c_set
      l[v] = k
      candidates.append(v)
    if len(candidates) == 0:
      recalculate_flag = False
      for v in g.vertices():
        if not v in cset:
          recalculate_flag = True
          print("fatal: vertex ",v," not in focused gflow ")
          # breakpoint()
      if recalculate_flag:
        return focus_gflow(g, gflow(g))
      for output in g.outputs(): #delete outputs from g function
        cset.pop(output)
      return [gf[0], cset, gf[2]] #[l, gflow, k]
    unprocessed.difference_update(set(candidates))
    k += 1

'''
Builds focused gflow graph as lists
'''
def build_line_rec(g, vertex, gflow):
  nlist = []
  for n in g.neighbors(vertex):
    if n in gflow[1] and vertex in gflow[1][n]:
      nlist.append(n)
  if not nlist:
    return [vertex]
  if len(nlist) > 1:
    # print("branching on vertex ",vertex)
    res = []
    for n in nlist:
      res.append(build_line_rec(g, n, gflow))
    return [vertex] + res
  else:
    return [vertex] + build_line_rec(g, nlist[0], gflow)

'''
Builds a Diagram out of the focused gflow graph for visualization purposes
'''
def build_diagram_from_lines(g, data):
  graph = GraphS()
  for node in range(0,len(data)):
    graph.add_vertex_indexed(data[node][0])
    recursive_diagram_build(graph, data[node][1:], data[node][0])
  graph.set_inputs(g.inputs())
  graph.set_outputs(g.outputs())
  for v in graph.vertices():
    graph.set_qubit(v, g.qubit(v))
    graph.set_row(v, g.row(v))
  return graph

'''
Recursive helper function for build_diagram_from_lines function 
'''
def recursive_diagram_build(graph, data, start):
  last_vertex = start
  already_visited = False
  for vertex in data:
    if isinstance(vertex, list):
      if not vertex[0] in graph.vertices():
        graph.add_vertex_indexed(vertex[0])
        recursive_diagram_build(graph, vertex[1:], vertex[0])
      # else:
      #   already_visited = True
      graph.add_edges([(last_vertex,vertex[0])])
      # Assumes this is the last vertex in data list
    else:
      if not vertex in graph.vertices():
        graph.add_vertex_indexed(vertex)
      else:
        already_visited = True
      graph.add_edges([(last_vertex,vertex)])
      last_vertex = vertex

    if already_visited:
      break

def get_odd_neighbourhood(g: BaseGraph[VT,ET], vertex_set):
  all_neighbors = set()
  for vertex in vertex_set:
    all_neighbors.update(set(g.neighbors(vertex)))
  odd_neighbors = []
  for neighbor in all_neighbors:
    if len(set(g.neighbors(neighbor)).intersection(vertex_set)) % 2 == 1:
      odd_neighbors.append(neighbor)
  return odd_neighbors


def gflow(g: BaseGraph[VT,ET]) -> Optional[Tuple[Dict[VT,int], Dict[VT,Set[VT]], int]]:
    l:     Dict[VT,int]      = {}
    gflow: Dict[VT, Set[VT]] = {}
    for v in g.outputs():
        l[v] = 0

    inputs = set(g.inputs())
    processed = set(g.outputs())
    vertices = set(g.vertices())
    k = 1
    while True:
        correct = set()
        #unprocessed = list()
        processed_prime = [v for v in processed.difference(inputs) if any(w not in processed for w in g.neighbors(v))]
        candidates = [v for v in vertices.difference(processed) if any(w in processed_prime for w in g.neighbors(v))]
        
        zerovec = Mat2([[0] for i in range(len(candidates))])
        #print(unprocessed, processed_prime, zerovec)
        m = bi_adj(g, processed_prime, candidates)
        cnot_maker = CNOTMaker()
        m.gauss(x=cnot_maker, full_reduce=True)
        for u in candidates:
            vu = zerovec.copy()
            vu.data[candidates.index(u)] = [1]
            x = get_gauss_solution(m, vu, cnot_maker.cnots)
            if x:
                correct.add(u)
                gflow[u] = {processed_prime[i] for i in range(x.rows()) if x.data[i][0]}
                l[u] = k

        if not correct:
            if not candidates:
                return l, gflow, k
            return None
        else:
            processed.update(correct)
            k += 1

def get_gauss_solution(gauss: Mat2, vec: Mat2, cnots: CNOTMaker):
  for cnot in cnots:
    vec.row_add(cnot.target,cnot.control)
  x = Mat2.zeros(gauss.cols(),1)
  for i,row in enumerate(gauss.data):
      got_pivot = False
      for j,v in enumerate(row):
          if v != 0:
              got_pivot = True
              x.data[j][0] = vec.data[i][0]
              break
      # zero LHS with non-zero RHS = no solutions
      if not got_pivot and vec.data[i][0] != 0:
          return None
  return x

#Copy from extract.py but necessary since otherwise we have circular import problems
def bi_adj(g: BaseGraph[VT,ET], vs:List[VT], ws:List[VT]) -> Mat2:
    """Construct a biadjacency matrix between the supplied list of vertices
    ``vs`` and ``ws``."""
    return Mat2([[1 if g.connected(v,w) else 0 for v in vs] for w in ws])
