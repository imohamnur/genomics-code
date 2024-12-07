from Bio import pairwise2
import matplotlib.pyplot as plt
import csv
import sys
def read_sequences_from_fastq(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        while True:
            header = file.readline().strip()
            sequence = file.readline().strip()
            plus = file.readline().strip()
            quality = file.readline().strip()
            if not sequence:
                break
            sequences.append(sequence)
    return sequences

def get_kmer_count_from_sequence(sequence, k=3, cyclic=False):
    kmers = {}
    for i in range(0, len(sequence)):
        kmer = sequence[i:i + k]
        length = len(kmer)
        if cyclic:
            if len(kmer) != k:
                kmer += sequence[:(k - length)]
        else:
            if len(kmer) != k:
                continue
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1
    return kmers

def get_debruijn_edges_from_kmers(kmers):
    edges = {}
    for k1 in kmers:
        for k2 in kmers:
            if k1 != k2:
                if k1[1:] == k2[:-1]:
                    edge = (k1[:-1], k2[:-1])
                    if edge in edges:
                        edges[edge] += 1  # Update edge multiplicity
                    else:
                        edges[edge] = 1
                if k1[:-1] == k2[1:]:
                    edge = (k2[:-1], k1[:-1])
                    if edge in edges:
                        edges[edge] += 1  # Update edge multiplicity
                    else:
                        edges[edge] = 1
    return edges

def trim_dead_ends(graph):
    # Find and remove nodes that have no out-edges (dead ends)
    dead_ends = [node for node in graph if not graph[node]]
    for node in dead_ends:
        del graph[node]
    return graph

def separate_contigs(graph):
    # Identify nodes with multiple out-edges (branching points)
    branching_points = [node for node in graph if len(graph[node]) > 1]
    contigs = []
    for point in branching_points:
        contig = []
        current = point
        while len(graph[current]) == 1:
            contig.append(current)
            next_node = graph[current][0]
            del graph[current]
            current = next_node
        contig.append(current)
        contigs.append(contig)
    return contigs

def find_eulerian_path(graph):
    def dfs(node):
        nonlocal circuit
        while graph[node]:
            neighbor = graph[node].pop()
            dfs(neighbor)
        circuit.append(node)

    start = list(graph.keys())[0]
    for node, edges in graph.items():
        if len(edges) > 0:
            start = node
            break

    circuit = []
    dfs(start)
    eulerian_path = circuit[::-1]
    return eulerian_path

def combine_kmers(sequences, k, cyclic=True):
    combined_kmers = {}
    for sequence in sequences:
        kmers = get_kmer_count_from_sequence(sequence, k, cyclic)
        for kmer, count in kmers.items():
            if kmer in combined_kmers:
                combined_kmers[kmer] += count
            else:
                combined_kmers[kmer] = count
    return combined_kmers

# Add a function for sequence alignment and plotting
def align_and_plot(assembled_genome, spike_protein_sequence):
    alignments = pairwise2.align.localms(assembled_genome, spike_protein_sequence, 2, -1, -0.5, -0.1)
    aligned_assembly = alignments[0][0]
    aligned_spike_protein = alignments[0][1]

    aligned_positions = []
    for i, (a, b) in enumerate(zip(aligned_assembly, aligned_spike_protein)):
        if a == b:
            aligned_positions.append((i, i))  # Store (assembly_coord, protein_coord)

    x = [pos[1] for pos in aligned_positions]
    y = [pos[0] for pos in aligned_positions]

    plt.scatter(x, y)
    plt.xlabel('Spike Protein Coordinates')
    plt.ylabel('Assembled Contig Coordinates')
    plt.title('Alignment of Assembled Contig to Spike Protein')
    plt.savefig('alignment_plot.png')  # Save as PNG file
    plt.show()


if __name__ == "__main__":
    file = 'output.fastq'
    sequences = read_sequences_from_fastq(file)

    all_kmers = combine_kmers(sequences, k=2, cyclic=False)
    print("K-mers:", all_kmers)

    edges = get_debruijn_edges_from_kmers(all_kmers)
    print("Edges:", edges)

    if edges:
        # Trimming dead-end nodes
        graph = {k: [] for k in set(sum(edges, ()))}
        for edge, count in edges.items():
            graph[edge[0]].append(edge[1])

        graph = trim_dead_ends(graph)
        print("Graph after trimming dead ends:", graph)

        # Separating contigs
        contigs = separate_contigs(graph)
        print("Contigs separated from branching points:", contigs)

        # Assembling genome from the Eulerian path
        eulerian_path = find_eulerian_path(graph)
        print("Eulerian Path:", eulerian_path)

        genome = ''.join([eulerian_path[i][0] for i in range(len(eulerian_path))])
        print("Assembled Genome:", genome)

        # Save the assembled genome to a file
        with open("assembled_genome.txt", "w") as output_file:
            output_file.write(genome)

        # Assuming you have the spike protein sequence stored in a variable
        spike_protein_sequence = "ACGT..."  # Replace with actual sequence

        # Align the assembled genome to the spike protein and plot
        align_and_plot(genome, spike_protein_sequence)
    else:
        print("No edges formed.")