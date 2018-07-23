import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from graphviz import Digraph
import pydot


class Vertex():
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.adj = defaultdict()


class Edge():
    def __init__(self, v1, v2):
        self.v1 = v1
        self.v2 = v2
        self.subst = list()


class Graph():
    def __init__(self):
        self.vertices = defaultdict()
        self.edges = defaultdict()

    def add_vertex(self, file):
        for record in SeqIO.parse(file, 'fasta'):
            self.vertices[(Vertex(record.id, str(record.seq)))] = []

    def adjacency_list(self):
        for vertex in self.vertices.keys():
            for vert in self.vertices.keys():
                if vert.name != vertex.name:
                    local_edge = Edge(vertex, vert)
                    for i in range(0, len(vertex.seq)):
                        if vertex.seq[i] != vert.seq[i]:
                            local_edge.subst.append(f"{vertex.seq[i]}{i+1}{vert.seq[i]}")
                    self.vertices[vertex].append({vert:local_edge})

    def visualize(self):
        vis = Digraph(comment='Cry toxin amino acid substitution graph')
        for i, j in self.vertices.items():
            vis.node(i.name, label=f"{i.name}")
            for num in j:
                for el in num.keys():
                    print(num[el].subst, '\n', i.name, el.name)
                    vis.edge(i.name, el.name, label=num[el].subst[0])
        vis.view()
        vis.save()


my_graph = Graph()
my_graph.add_vertex('/home/reverend_casy/toxin_dummy_aligned.fasta')
my_graph.adjacency_list()

if list(my_graph.vertices[list(my_graph.vertices.keys())[0]][0].keys())[0].name == list(my_graph.vertices.keys())[1].name:
    print("excessive edge!")




my_graph.visualize()
