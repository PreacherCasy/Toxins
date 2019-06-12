from Bio import AlignIO
import numpy as np
from Bio import SeqRecord
from collections import defaultdict
from graphviz import Graph
from collections import deque
from itertools import chain
import pandas as pd
import subprocess
import argparse


def merge_indels(label):
    dum_list = list()
    i = 0
    while i < len(label):
        if i < len(label) - 1 and int(label[i][1]) == int(label[i][1]) + 1:
            if i < len(label) - 1 and label[i][2] == '-' and label[i+1][0] == '-':
                dum_list.append(f"{label[i][0:2]}|{label[i+1][1:]}")
                i += 2
            elif i < len(label) - 1 and label[i][0] == '-' and label[i+1][2] == '-':
                dum_list.append(f"{label[i][2]}{label[i][1]}|{label[i+1][1]}{label[i+1][0]}")
                i += 2
            else:
                dum_list.append(label[i])
                i += 1
        else:
            dum_list.append(label[i])
            i += 1
    return dum_list


class Vertex:

    def __init__(self, vert):
        self.name = vert.id
        self.seq = str(vert.seq)
        self.in_edges = dict()
        self.out_edges = dict()

    def __str__(self):
        return self.name


class Edge:

    def __init__(self, k1, k2):
        self.start = k1
        self.finish = k2
        self.label = list()
        self.status = 1
        self.coords = (0, 0)

    def elongate_label(self, subst):
        self.label.append(subst)

    def merge_label(self):
        self.label = merge_indels(self.label)

    def __str__(self):
        joined_label = ','.join(self.label)
        return f"{self.start.name}->{self.finish.name}, edge label: {joined_label}"


class AAGraph:
    def __init__(self, alignment):
        self.alignment = open(alignment)
        self.al = AlignIO.read(self.alignment, 'fasta')
        self.adjMat = np.triu(np.full((len(self.al), len(self.al)), 1))
        np.fill_diagonal(self.adjMat, 0)
        self.distMat = np.zeros((len(self.al), len(self.al)))
        self.alignment.close()
        self.vertices = list()
        self.edges = list()
        self.filtered_edges = list()

    def upload_vertices(self, start, stop, entry_name):
        if not start:
            start = 0
         if not stop:
            stop = len(self.al[0].seq)
        indices = []
        for entry in self.al:
            if entry.id == entry_name:
                query = str(entry.seq).replace('-', '')[start:stop]
                for i in range(1, len(str(entry.seq))):
                    if str(entry.seq)[i] != '-':
                        for j in range(i, len(str(entry.seq))):
                            if str(entry.seq)[i:j].replace('-', '') == query:
                                indices = [i, j]
        for entry in self.al:
            self.vertices.append(Vertex(SeqRecord.SeqRecord(id=entry.id, seq=str(entry.seq)[indices[0]:indices[1]])))

    def matBuild(self):
        self.distMat = np.zeros((len(self.al), len(self.al)))
        for k in range(1, len(self.vertices)):
            for i in range(k):
                local_edge = Edge(self.vertices[k], self.vertices[i])
                for j in range(0, len(self.vertices[i].seq)):
                    if self.vertices[k].seq[j] != self.vertices[i].seq[j]:
                        self.distMat[i][k] += 1
                        local_edge.elongate_label(f"{self.vertices[k].seq[j]}{str(j + 1)}{self.vertices[i].seq[j]}")
                local_edge.coords = (i, k)
                local_edge.merge_label()
                self.edges.append(local_edge)
                self.vertices[k].out_edges[self.vertices[i]] = local_edge
                self.vertices[i].in_edges[self.vertices[k]] = local_edge

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
            for edge in self.edges:
                if self.adjMat[edge.coords[0], edge.coords[1]] == 0:
                    edge.status = 0

    def non_gapped(self, root, threshold):

        trans_dict = defaultdict(list)
        final_dict = dict()
        for edge in self.edges:
            trans_dict[edge.start].append((edge.finish, edge.label))
            trans_dict[edge.finish].append((edge.start, list(map(lambda x: x[::-1], edge.label))))

        visited, queue = set(), deque([x for x in self.vertices if x.name == root])
        while queue:
            vertex = queue.popleft()
            save_list = list()
            for neighbour in trans_dict[vertex]:
                if neighbour[0] not in visited:
                    visited.add(neighbour[0])
                    if (''.join(neighbour[1]).count('-') == -1 or ''.join(neighbour[1]).count('-') <= threshold) \
                            and neighbour not in save_list:
                        save_list.append(neighbour)
            if len(save_list) > 0:
                if vertex not in final_dict.keys():
                    final_dict[vertex] = save_list
                else:
                    final_dict[vertex] = list(chain(final_dict[vertex], save_list))
            for neighbour in save_list:
                queue.append(neighbour[0])
        new_edges = list()
        for item in final_dict.items():
            for finish in item[1]:
                new_edge = Edge(item[0], finish[0])
                new_edge.label = finish[1]
                new_edges.append(new_edge)
        self.filtered_edges = new_edges

    def export(self, title):
        export = pd.DataFrame({'Start_vertex': list(map(lambda x: x.start, self.edges)), \
                               'End_vertex': list(map(lambda x: x.finish, self.edges)), \
                               'Mismatches': list(map(lambda x: x.label, self.edges))})
        export.to_csv(title, sep='\t', header=True)
        return export

    def visualize(self, *args):
        vis = Graph(comment='Cry toxin amino acid substitution graph', strict=True)
        build_graph = self.filtered_edges if len(self.filtered_edges) > 0 else self.edges

        for edge in build_graph:
            if edge.status == 1:
                vis.node(edge.start.name, label=f"{edge.start.name}")
                vis.node(edge.finish.name, label=f"{edge.finish.name}")
                vis.edge(edge.start.name, edge.finish.name, label=", ".join(edge.label))
        if args:
            vis.render(args[0], view=True)
        vis.view()
        vis.save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='graph_visualization')
    parser.add_argument('-ai', help='Name/path to your fasta file', metavar='Alignment input',
                        type=str, required=True)
    parser.add_argument('-ao', help='Name/path to your alignment file', metavar='Alignment output',
                        type=str, required=True)
    parser.add_argument('-n', help='Name of the defining sequence', metavar='Definin sequence',
                        type=str, required=True)
    parser.add_argument('-s', help='Starting index of substring to be converted into vertex', metavar='Starting index',
                        type=int, required=True)
    parser.add_argument('-f', help='Ending index of substring to be converted into vertex', metavar='Ending index',
                        type=int, required=True)
    parser.add_argument('-ng', help='Defines the non-gapped subgraph mode', required=False, action='store_true')
    parser.add_argument('-r', help='Root vertex for non-gapped graph reconstruction', metavar='Non-gapped mode root',
                        type=str, required=False)
    parser.add_argument('-thr', help='Maximum number of gaps allowed to maintain an edge in non-gapped graph',
                        metavar='Non-gapped mode indel threshold', nargs='?', type=int, required=False)
    parser.add_argument('-csv', help='Path to the csv file', metavar='csv path',
                        type=str, required=False)
    parser.add_argument('-gv', help='Path to the dot file', metavar='dot path',
                        type=str, required=False)
    parser.set_defaults(feature=True)

    args = parser.parse_args()
    ai, ao, n, s, f, ng, r, thr, csv, gv = args.ai, args.ao, args.n, args.s, \
                                           args.f, args.ng, args.r, args.thr, args.csv, args.gv

    subprocess.call(f"mafft --quiet {ai} > {ao}", shell=True)
    g = AAGraph(ao)
    g.upload_vertices(s, f, n)
    g.matBuild()
    g.TR()
    if ng:
        if thr:
            g.non_gapped(r, thr)
        else:
            g.non_gapped(r)
    if csv:
        g.export(csv)
    if gv:
        g.visualize(gv)
    else:
        g.visualize()
