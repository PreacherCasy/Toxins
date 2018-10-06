from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
import numpy as np
from scipy.sparse.csgraph import connected_components as cc
from graphviz import Graph
import argparse
import subprocess

class graph():
    def __init__(self, alignment):
        self.alignment = alignment
        self.distmat = np.triu(DistanceCalculator('identity').get_distance(AlignIO.read(open(alignment), 'fasta')),1)
        self.adj_matrix = np.triu(np.full((len(AlignIO.read(open(alignment), 'fasta')), len(AlignIO.read(open(alignment),
                                                                                                'fasta'))), 1), 1)
        self.vertices = dict()

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
    parser.set_defaults(feature=True)

    args = parser.parse_args()
    ai, ao = args.ai, args.ao

    subprocess.call(f"mafft {ai} > {ao}", shell=True)
    g = graph(ao)
    g.edge_cutter()
    g.visualize()
