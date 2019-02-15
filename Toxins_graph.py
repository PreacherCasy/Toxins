from Bio import AlignIO
import numpy as np
from collections import defaultdict
from graphviz import Graph
import argparse
import subprocess


class AAGraph():
    def __init__(self, alignment):
        self.alignment = open(alignment)
        self.al = AlignIO.read(self.alignment, 'fasta')
        self.adjMat = np.triu(np.full((len(self.al), len(self.al)), 1))
        np.fill_diagonal(self.adjMat, 0)
        self.distMat = np.zeros((len(self.al), len(self.al)))
        self.alignment.close()
        self.vertices = list()
        self.edges = defaultdict(list)

    def upload_vertices(self):
        for entry in self.al:
            self.vertices.append(entry.id)

    def matBuild(self):
        self.distMat = np.zeros((len(self.al), len(self.al)))
        for k in range(1, len(self.al)):
            for i in range(k):
                for j in range(0, len(self.al[i].seq)):
                    if self.al[k].seq[j] != self.al[i].seq[j]:
                        self.distMat[i][k] += 1
                        self.edges[(self.al[k].id, self.al[i].id)].\
                            append(f"{self.al[k].seq[j]}{str(j + 1)}{self.al[i].seq[j]}")
        self.distMat = self.distMat

    def TR(self):
        for k in range(1, len(self.distMat)):
            for i in range(k):
                for j in range(k + 1, len(self.distMat)):
                    viaK = max(self.distMat[i, k], self.distMat[k, j])
                    if viaK < abs(self.distMat[i, j]) and abs(self.distMat[i, j])!= 0:
                        self.distMat[i, j] = -viaK
            for i in range(len(self.distMat)):
                for j in range(len(self.distMat)):
                    if self.distMat[i, j] < 0:
                        self.adjMat[i, j] = 0

    def visualize(self):
        vis = Graph(comment='Cry toxin amino acid substitution graph', strict=True)
        for start, i in enumerate(self.vertices):
            vis.node(i, label=f"{i}")
            for end in range(start, ):
                if self.adjMat[end][start] == 1:
                    item = (i, self.vertices[end])
                    vis.edge(i, self.vertices[end], label=", ".join(self.edges[item]))
        vis.view()
        vis.save()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='graph_viualization')
    parser.add_argument('-ai', help='Name/path of your fasta file', metavar='File',
                        type=str, required=True)
    parser.add_argument('-alg', help='Preferred aligner name and/or absolute path', metavar='Program',
                        type=str, required=True)
    parser.add_argument('-ao', help='Name/path to your alignment file', metavar='File',
                        type=str, required=True)
    parser.set_defaults(feature=True)

    args = parser.parse_args()
    ai, alg, ao = args.ai, args.alg, args.ao
    subprocess.run(f"{alg} --quiet {ai} > {ao}", shell=True)
    g = AAGraph(ao)
    g.upload_vertices()
    g.matBuild()
    g.TR()
    g.visualize()
