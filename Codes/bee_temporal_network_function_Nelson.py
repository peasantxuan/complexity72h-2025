from math import log2

import matplotlib.patches as mpatches
# import matplotlib.pyplot as plt
import networkx as nx
import networkx_temporal as tx
import numpy as np
import pandas as pd

# https://github.com/hfelippe/network-MI
from hfelippe.network_MI.functions import NMI, DCNMI


tab10 = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]


def load_labels(filename, folder='.'):
    """ Load labels from a CSV file, skipping the header. """
    with open(f"{folder}/{filename}", 'r') as f:
        return [x.strip('"') for x in f.readline().strip().split(',')[1:]]


def load_vector(filename, folder='.', quotechar=None, usecols=None, dtype=float):
    """ Load a vector from a CSV file, skipping the header. """
    with open(f"{folder}/{filename}", 'r') as f:
        cols = len(f.readline().strip().split(','))
    if cols == 1:
        return np.loadtxt(f"{folder}/{filename}", skiprows=1)
    return np.loadtxt(f"{folder}/{filename}", delimiter=',', skiprows=1, dtype=dtype,
                      quotechar=quotechar, usecols=usecols if usecols is not None else range(1, cols))


def load_matrix(filename, folder='.'):
    """ Load a matrix from a CSV file, skipping the first row and column. """
    with open(f"{folder}/{filename}", 'r') as f:
        cols = len(f.readline().strip().split(','))
    return np.loadtxt(f"{folder}/{filename}", delimiter=',', skiprows=1,
                      usecols=range(1, cols))


def build_matrices(A, B, P):
    """ Build temporal matrices from base matrix A, change matrix B, and environmental factors P. """
    matrices = []
    for i in range(len(P)):
        matrix = A + (B * P[i])
        matrices.append(matrix)
    return matrices


def apply_threshold(matrix):
    """ Apply a threshold to the matrix to filter out small values. """
    # values = matrix.flatten()
    offdiag = matrix[np.where(~np.eye(matrix.shape[0],dtype=bool))]
    mu = np.mean(offdiag)
    sigma = np.std(offdiag)
    pos = np.where(matrix >= mu + sigma / 2, matrix, 0)
    neg = np.where(matrix <= mu - sigma / 2, matrix, 0)
    mat = pos + neg
    np.fill_diagonal(mat, 0)  # Remove self-loops
    return mat


def scale_width(weights):
    """Scale edge widths based on their weights. """
    return [
        max(1.2, np.abs(w * 6.8))
        for w in weights
    ]


# def color_signed_edge(weight, blue="#1f77b4", red="#d62728"):
def color_signed_edge(weight, blue="#1f77b4", red="#d62728"):
    """ Color edges based on their weight. """
    if weight > 0:
        return blue
    if weight < 0:
        return red
    return 'gray'


def build_temporal_graph(matrices):
    """ Build a temporal graph from a list of matrices. """
    TG = tx.from_snapshots([
        nx.from_numpy_array(apply_threshold(matrices[i]), create_using=nx.DiGraph())
        for i in range(len(matrices))
    ])

    # Remove self-loops.
    for G in TG:
        G.remove_edges_from([(u, u) for u in G.nodes() if G.has_edge(u, u)])

    return TG


def draw_graph(TG, names=None, labels=[], title=''):
    """ Draw a temporal graph from a list of matrices."""
    # Compute node positions.
    pos = nx.circular_layout(TG[0])

    # Number of rows and columns for subplots.
    nrows = 2
    while True:
        if len(TG) % nrows == 0:
            ncols = len(TG) // nrows
            break
        nrows += 1

    # Node and edge colors/sizes.
    # node_color = [
    #     tab10[i % len(tab10)]
    #     for i in range(len(labels))
    # ]
    node_color = "#777F92"
    edge_color = [
        [color_signed_edge(w) for u, v, w in G.edges(data="weight")]
        for G in TG
    ]
    width = [
        scale_width([w for u, v, w in G.edges(data="weight")])
        for G in TG
    ]

    # Draw graphs.
    fig = tx.draw(
        TG,
        pos=pos,
        nrows=nrows,
        ncols=ncols,
        figsize=(10, 5),
        node_size=500,
        node_color=node_color,
        temporal_edge_color=edge_color,
        temporal_width=width,
        with_labels=False,
        labels={i: labels[i] for i in range(len(labels))},
        connectionstyle="arc3,rad=0.1",
        # suptitle="BEEFUN",
        names=names,
    )

    # Prepare legend (option 1).
    # line1 = Line2D([], [], color="white", marker='o', markersize=10, markerfacecolor="slategray")
    # ...
    # plt.legend((line, line, line, line), ('Thing 1', 'Thing 2', 'Thing 3', 'Thing 4'), numpoints=1, loc=1)

    # Prepare legend (option 2).
    # handles = [
    #     mpatches.Patch(color=tab10[i % len(tab10)], label=labels[i])
    #     for i in range(len(labels))
    # ]

    # Add legend for species.
    # fig.legend(
    #     handles=handles,
    #     labels=labels,
    #     loc='lower center',
    #     bbox_to_anchor=(0.5, 1.05),
    #     ncol=5,
    #     fontsize='small',
    #     title=title or 'Species',
    #     title_fontsize='medium',
    #     frameon=False,
    # )

    return fig


def graph_entropy(G):
    m = G.size()
    deg = dict(nx.degree(G))
    H = sum(deg[v]/m * log2(deg[v]/m) for v in G.nodes)
    return H


def compute_centralizations(TG):
    """ Calculate centralization for degree, in_degree, and out_degree. """

    centralizations = {}

    for degree in ("degree", "in_degree", "out_degree"):
        centralizations[degree] = []

        for G in tx.from_multigraph(TG):
            centrality = getattr(G, degree)()

            maximum = G.order() - 1  # Highest possible degree.
            minimum = 1              # Minimum possible degree.

            # Highest theoretical sum of values in a simple star-like graph.
            scalar = sum(maximum - minimum for n in range(G.order()-1))

            # List of node degree centrality values.
            centrality = [d for n, d in centrality]
            max_centrality = max(centrality)
            centralization = sum(max_centrality - c for c in centrality)

            # Normalize the centralization value.
            centralization /= (scalar or 1)

            centralizations[degree].append(centralization)

    return centralizations


def compute_biedges(TG, years):
    """ Compute the number of bidirectional edges in a temporal graph. """
    # Calculate bidirectional edge ratio.
    biedges = []

    for i, G in enumerate(tx.from_multigraph(TG)):
        biedge = sum([1 if G.has_edge(u, v) and G.has_edge(v, u) else 0 for u in G.nodes() for v in G.nodes() if u != v])
        biedge_pos = sum([1 if G.has_edge(u, v) and G.has_edge(v, u) and G[u][v]["weight"] > 0 else 0 for u in G.nodes() for v in G.nodes() if u != v])
        biedge_neg = sum([1 if G.has_edge(u, v) and G.has_edge(v, u) and G[u][v]["weight"] < 0 else 0 for u in G.nodes() for v in G.nodes() if u != v])
        theoretical_max = (G.order() * (G.order() - 1)) / 2  # Maximum number of edges in a complete graph.
        print(f"Bidirectional edges ({years[i]}): {biedge} ({biedge/theoretical_max:.2%}), positive: {biedge_pos}, negative: {biedge_neg}")
        biedges.append(biedge / theoretical_max)

    return biedges


def compute_nmis(TG, dc=False):
    """ Compute Normalized Mutual Information (NMI)
    between consecutive years in a temporal graph,
    optionally employing degree correction (dc). """
    nmis = []

    for i in range(len(TG)):
        if i > 0:
            n = TG[i].order()
            e1 = pd.Index(list(TG[i].edges()))
            e2 = pd.Index(list(TG[i-1].edges()))
            nmi = (DCNMI if dc else NMI)(n, e1, e2)
            print(f"{'DC' if dc else ''}NMI between year {2015+i-1} and {2015+i}: {nmi}")
            nmis.append(nmi)

    return nmis
