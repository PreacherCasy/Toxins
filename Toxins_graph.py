from Bio import SeqIO
from graphviz import Digraph
import argparse

class Vertex():
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
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
        self.adj_list = dict()

    def add_vertex(self, file):
        for record in SeqIO.parse(file, 'fasta'):
            self.vertices.append(Vertex(record.id, str(record.seq)))


    def find_path(self, start, end, path):
        try:
            path = path + [start]
            if start == end:
                return True
            for node in self.adj_list[start]:
                if node not in path:
                    newpath = self.find_path(node, end, path)
                    if newpath: return True
            return False
        except:
            return False

    def checker(self, start, end, path):
        if self.find_path(start, end, path) == True:
            return False
        else:
            return True



    def adjacency_list(self):
        for vertex in self.vertices:
            self.adj_list[vertex] = list()
            for vert in self.vertices:
                if vert.name != vertex.name and self.checker(vert, vertex, path = []) is True and vertex not in vert.out_edges.keys():
                    local_edge = Edge(vertex, vert)
                    for i in range(0, len(vertex.seq)):
                        if vertex.seq[i] != vert.seq[i]:
                            local_edge.subst.append(f"{vertex.seq[i]}{i+1}{vert.seq[i]}")
                    if len(local_edge.subst) < vertex.baseline:
                        vertex.out_edges.clear()
                        vertex.out_edges[vert], vertex.baseline = local_edge, len(local_edge.subst)
                        self.adj_list[vertex].clear()
                        self.adj_list[vertex].append(vert)
                    elif len(local_edge.subst) == vertex.baseline or vertex.baseline == 0:
                        vertex.out_edges[vert], vertex.baseline = local_edge, len(local_edge.subst)
                        self.adj_list[vertex].append(vert)
        print(self.adj_list)

    def visualize(self):
        vis = Digraph(comment='Cry toxin amino acid substitution graph', strict=True)
        for i in self.vertices:
            vis.node(i.name, label=f"{i.name}")
            for item in i.out_edges.items():
                vis.edge(i.name, item[0].name, label = ", ".join(item[1].subst))

        vis.view()
        vis.save()

my_graph = graph()
my_graph.add_vertex('/home/reverend_casy/toy_set.fasta')
my_graph.adjacency_list()
my_graph.visualize()
