import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from graphviz import Digraph
import argparse



class Vertex():
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.adj = dict()

class Edge():
    def __init__(self, v1, v2):
        self.v1 = v1
        self.v2 = v2
        self.subst = list()

class graph():
    def __init__(self):
        self.vertices = dict()
        self.edges = dict()

    def add_vertex(self, file):
        for record in SeqIO.parse(file, 'fasta'):
            self.vertices[(Vertex(record.id, str(record.seq)))] = []

    def adjacency_list(self):
        for vertex in self.vertices.keys():

            for vert in list(self.vertices.keys())[list(self.vertices.keys()).index(vertex):]:
                print(vert.name)
                if vert.name != vertex.name:
                    local_edge = Edge(vertex, vert)
                    for i in range(0, len(vertex.seq)):
                        if vertex.seq[i] != vert.seq[i]:
                            local_edge.subst.append(f"{vertex.seq[i]}{i+1}{vert.seq[i]}")
#                    if len(local_edge.subst) <= k:
                    try:
                        if len(local_edge.subst) < len(self.vertices[vertex][0][1].subst):
                            self.vertices[vertex][0] = [vert, local_edge]
                        elif len(local_edge.subst) == len(self.vertices[vertex][0][1].subst):
                            self.vertices[vertex].append([vert, local_edge])
                    except:
                        self.vertices[vertex].append( [vert, local_edge])


    def component(self):
        vis = Digraph(comment='Cry toxin amino acid substitution graph')
        for i, j in self.vertices.items():
            vis.node(i.name, label=f"{i.name}")
            for item in j:
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
