from Bio import SeqIO
from Bio import AlignIO
import numpy as np
from collections import defaultdict
from graphviz import Graph
import argparse
import subprocess


class AAGraph():
    def __init__(self, alignment):
        self.alignment = alignment
        self.al = AlignIO.read(open(self.alignment), 'fasta')
        self.adjMat = np.triu(np.full((len(AlignIO.read(open(alignment), 'fasta')), len(AlignIO.read(open(alignment), 'fasta'))), 1))
        np.fill_diagonal(self.adjMat, 0)
        self.distMat = np.zeros((len(AlignIO.read(open(alignment), 'fasta')), len(AlignIO.read(open(alignment), 'fasta'))))
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

dum = AAGraph('/home/reverend_casy/Bt_toxins/toy_set_aligned.fasta')
dum.upload_vertices()
dum.matBuild()
dum.TR()
dum.visualize()
