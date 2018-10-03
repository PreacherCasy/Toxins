from Bio import SeqIO
import numpy as np
from scipy.sparse.csgraph import connected_components as cc
from graphviz import Graph
import argparse
import subprocess

class graph():
    def __init__(self, distances, alignment):
        self.distances = distances
        self.alignment = alignment
        self.adj_matrix = np.triu(np.full((11,11), 1), 1)
        self.vertices = dict()
        self.edges = list()
        self.distmat = np.zeros((11, 11))


    def populate_distances(self):
        filling_list = []
        with open(self.distances) as handle:
            for line in handle:
                line = list(map(float, filter(None, line.replace('\t', ' ').replace('\n', ' ').split(' '))))
                filling_list.append(line)
        for i in range(0, len(filling_list)):
            for j in range(0, len(filling_list[i])):
                self.distmat[i, -j] = filling_list[i][-j]

    def edge_cutter(self):
        if np.count_nonzero(self.distmat) > 0:
            i, j = np.argwhere(self.distmat == np.amax(self.distmat))[0][0], np.argwhere(self.distmat
                                                                                             == np.amax(self.distmat))[0][1]
            self.adj_matrix[i][j] = 0
            if cc(self.adj_matrix)[0] != 1:
                self.adj_matrix[i][j] = 1
            self.distmat[i][j] = 0
            self.edge_cutter()

    def upload_vertices(self):
        vert_list, self.g = list(), dict()
        for record in SeqIO.parse(self.alignment, 'fasta'):
            vert_list.append(record)

        for i in range(0, len(vert_list)):
            self.g[vert_list[i].id] = []
            for j in range(0, len(self.adj_matrix[i, :])):
                label = list()
                if self.adj_matrix[i, j] == 1:
                    for pos in range(0, len(vert_list[j].seq)):
                        if vert_list[i].seq[pos] != vert_list[j].seq[pos]:
                            label.append(f"{vert_list[i].seq[pos]}{pos + 1}{vert_list[j].seq[pos]}")
                    self.g[vert_list[i].id].append((vert_list[j].id, label))

    def visualize(self):
        self.upload_vertices()
        vis = Graph(comment='Cry toxin amino acid substitution graph', strict=True)
        for i in self.g.keys():
            vis.node(i, label=f"{i}")
            for item in self.g[i]:
                vis.edge(i, item[0], label=", ".join(item[1]))

        vis.view()
        vis.save()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='graph_viualization')
    parser.add_argument('-ai', help='Name/path to your fasta file', metavar='File',
                        type=str, required=True)
    parser.add_argument('-ao', help='Name/path to your alignment file', metavar='File',
                        type=str, required=True)
    parser.add_argument('-d', help='Name/path to your distance matrix file', metavar='File',
                        type=str, required=True)
    # parser.add_argument('-ms', help='MAFFT settings',
    #                     type=str, required=False, default=None)
    parser.set_defaults(feature=True)

    args = parser.parse_args()
    d, ai, ao = args.d, args.ai, args.ao

    subprocess.call(f"mafft {ai} > {ao}", shell=True)
    g = graph(d, ao)
    g.populate_distances()
    g.edge_cutter()
    g.visualize()
