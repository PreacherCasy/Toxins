from Bio import SeqIO
from collections import defaultdict
from graphviz import Digraph
import argparse

class Vertex():
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.in_edges = dict()
        self.out_edges = dict()
        self.baseline = 0

class Edge():
    def __init__(self, v1, v2):
        self.v1 = v1
        self.v2 = v2
        self.subst = list()

class graph():
    def __init__(self):
        self.vertices = list()
        self.outer_set = set()
        self.edges = set()

    def add_vertex(self, file):
        for record in SeqIO.parse(file, 'fasta'):
            self.vertices.append(Vertex(record.id, str(record.seq)))

    def acyclic(self, node1, node2):
        self.outer_set.add(node2)
        if node1 in node2.out_edges:
            return False
        else:
            for i in node2.out_edges.keys():
                if i not in self.outer_set:
                    self.acyclic(node1, i)

    def adjacency_list(self):
        for vertex in self.vertices:
            for vert in self.vertices:
                if vert.name != vertex.name  and self.acyclic(vert, vertex) != False and vertex not in vert.out_edges.keys():
                    local_edge = Edge(vertex, vert)
                    for i in range(0, len(vertex.seq)):
                        if vertex.seq[i] != vert.seq[i]:
                            local_edge.subst.append(f"{vertex.seq[i]}{i+1}{vert.seq[i]}")
                    if len(local_edge.subst) < vertex.baseline:
                        vertex.out_edges.clear()
                        vertex.out_edges[vert] = local_edge
                        vertex.baseline = len(local_edge.subst)
                        vert.in_edges.clear()
                        vert.in_edges[vertex] = local_edge
                    elif len(local_edge.subst) == vertex.baseline or vertex.baseline == 0:
                        vertex.out_edges[vert] = local_edge
                        vertex.baseline = len(local_edge.subst)
                        vert.in_edges[vertex] = local_edge

    def visualize(self):
        vis = Digraph(comment='Cry toxin amino acid substitution graph', strict=True)
        for i in self.vertices:
            vis.node(i.name, label=f"{i.name}")
            for item in i.out_edges.items():
                vis.edge(i.name, item[0].name, label = ", ".join(item[1].subst))

        vis.view()
        vis.save()
        
if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Cry toxin amino acid substitution graph')
        parser.add_argument('-k', help='Maximal number of mismatches for adjacency', metavar='Int', type=int, default=1)
        parser.add_argument('-p', help='add path to a file', metavar='Str', type=str)
        args = parser.parse_args()

        p = args.p

        my_graph = graph()
        my_graph.add_vertex(f"{p}")
        my_graph.adjacency_list()
        my_graph.component()
