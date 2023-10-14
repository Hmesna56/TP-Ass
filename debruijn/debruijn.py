#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Hamza Mesnaoui"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Hamza Mesnaoui"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Hamza Mesnaoui"
__email__ = "hamzamesnaoui1@gmail.com"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    with open(fastq_file, 'r') as file:
        for line in file:
            yield next(file).strip()
            next(file)
            next(file)


def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    read_length = len(read)
    for i in range(read_length - kmer_size + 1):
        kmer = read[i:i + kmer_size]
        #print(kmer)
        yield kmer


def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}

    read_sequences = read_fastq(fastq_file)

    for sequence in read_sequences:

        kmer_generator = cut_kmer(sequence, kmer_size)
    
        for kmer in kmer_generator:
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    print(kmer_dict)
    return kmer_dict


def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    graph = nx.DiGraph()
    for kmer, occurrence in kmer_dict.items():
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph.add_edge(prefix, suffix, weight=occurrence)
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    #modified_graph = graph.copy()
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)
        elif delete_entry_node:
            graph.remove_nodes_from(path[:-1])
            # Remove the first node
        elif delete_sink_node:
            graph.remove_nodes_from(path[1:])
            # Remove the last node
        else:
            graph.remove_nodes_from(path[1:-1])
    return graph


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    '''weight_stdev = statistics.stdev(weight_avg_list)
    if weight_stdev > 0:
        best_path_index = weight_avg_list.index(max(weight_avg_list))
        del path_list[best_path_index]
    else:
        length_stdev = statistics.stdev(path_length)
        if length_stdev > 0:
            longest_path_index = path_length.index(max(path_length))
            del path_list[longest_path_index]
        else:
            best_path_index = random.randint(0, len(path_list) - 1)
            del path_list[best_path_index]

    return remove_paths(graph, path_list, delete_entry_node, delete_sink_node)'''

    if not path_list:
        return graph  # Aucun chemin à sélectionner, retournez le graphe tel quel

    if len(path_list) == 1:
        return graph  # Il n'y a qu'un seul chemin, retournez le graphe tel quel

    std_weight = statistics.stdev(weight_avg_list)
    std_length = statistics.stdev(path_length)

    if std_weight > 0:
        # Sélectionnez le chemin avec le poids moyen le plus élevé
        num_best_path = weight_avg_list.index(max(weight_avg_list))
    elif std_length > 0:
        # Sélectionnez le chemin avec la longueur maximale
        num_best_path = path_length.index(max(path_length))
    else:
        # Choisissez aléatoirement un chemin
        num_best_path = random.randint(0, len(path_list) - 1)

    del path_list[num_best_path]
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph

def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    all_paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))

    path_lengths = [len(path) for path in all_paths]
    path_weights = [path_average_weight(graph, path) for path in all_paths]

    graph = select_best_path(graph, all_paths, path_lengths, path_weights, 
                            delete_entry_node=False, delete_sink_node=False)

    return graph

def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble = False

    for node in graph.nodes():
        list_predecesseur = list(graph.predecessors(node))

        if len(list_predecesseur) > 1:
            for i in range(len(list_predecesseur) - 1):
                for j in range(i + 1, len(list_predecesseur)):
                    ancestor_node = nx.lowest_common_ancestor(graph, list_predecesseur[i], list_predecesseur[j])

                    if ancestor_node is not None:
                        bubble = True
                        break  # Exit the inner loop

                if bubble:
                    break  # Exit the outer loop

    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, ancestor_node, node))

    return graph

def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    entry_tips = []
    for node in graph.nodes():
        if node in starting_nodes:
            continue  # Skip starting nodes

        predecessors = list(graph.predecessors(node))
        if len(predecessors) <= 1:
            continue  # Skip nodes with one or zero predecessors

        for start_node in starting_nodes:
            if nx.has_path(graph, start_node, node):
                entry_tips.extend(nx.all_simple_paths(graph, start_node, node))

    if len(entry_tips) > 1:
        list_weights = [path_average_weight(graph, path) for path in entry_tips]
        list_lengths = [len(path) for path in entry_tips]
        best_path = select_best_path(graph, entry_tips, list_lengths, list_weights, delete_entry_node=True, delete_sink_node=False)
        graph = solve_entry_tips(best_path, starting_nodes)

    return graph

def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    exit_tips = []

    for node in graph:
        successors = list(graph.successors(node))

        if len(successors) <= 1 or node in ending_nodes:
            continue

        for end_node in ending_nodes:
            if nx.has_path(graph, source=node, target=end_node):
                exit_tips.extend(nx.all_simple_paths(graph, source=node, target=end_node))

        if exit_tips:
            break

    if len(exit_tips) > 1:
        list_lengths = [len(path) for path in exit_tips]
        list_weights = [path_average_weight(graph, path) for path in exit_tips]
        best_path = select_best_path(graph, exit_tips, list_lengths, list_weights, delete_entry_node=False, delete_sink_node=True)
        graph = solve_out_tips(best_path, ending_nodes)

    return graph

def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    return [node for node in graph if not any(graph.predecessors(node))]

def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    return [node for node in graph if not any(graph.successors(node))]

def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs = []
    
    for start_node in starting_nodes:
        for end_node in ending_nodes:
            if nx.has_path(graph, start_node, end_node):
                for path in nx.all_simple_paths(graph, start_node, end_node):
                    contig = path[0]  # Le premier nœud d'entrée
                    for node in path[1:]:
                        contig += node[-1]  # Ajouter le dernier nucléotide du nœud
                    contig_length = len(contig)
                    contigs.append((contig, contig_length))
    
    return contigs

   
def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """

    with open(output_file, 'w') as file:  # Open the file in write mode ('w')
        for i, (sequence, length) in enumerate(contigs_list):
            # Create a FASTA entry for the current contig
            header = f">Contig_{i + 1} Length={length}\n"
            file.write(header)
            file.write(sequence + '\n')

def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    '''for a in  read_fastq(args.fastq_file):
        print(list(cut_kmer(a, 4)))'''
    dict_reads = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(dict_reads)
    graph = simplify_bubbles(graph)
    graph = solve_entry_tips(graph, get_starting_nodes(graph))
    graph = solve_out_tips(graph, get_sink_nodes(graph))
    #remove_paths(build_graph(build_kmer_dict(args.fastq_file, 4)), path_list, delete_entry_node, delete_sink_node)
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Resolve bubbles and get the paths to remove
    
    L_contigs = get_contigs(graph, get_starting_nodes(graph), get_sink_nodes(graph))
    save_contigs(L_contigs, args.output_file)



if __name__ == '__main__': # pragma: no cover
    main()
