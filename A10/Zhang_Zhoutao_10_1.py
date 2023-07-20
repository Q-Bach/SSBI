import argparse
import numpy as np
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import powerlaw
import sys
import warnings
warnings.filterwarnings("ignore")

def DFS(graph, node, temp, visited):
    # Helper function to perform depth first search
    visited[node] = 1
    temp = np.append(temp, node)
    for i in graph[node]:
        if (visited.get(i) is None):
            temp, visited = DFS(graph, i, temp, visited)
    return temp, visited

def main():
    # set up an argument parser 
    args_parser = argparse.ArgumentParser(
        prog="Zhang_Zhoutao_10_01",
        description="Protein Protein Interreaction Network")
    args_parser.add_argument(
        "-i",
        "--input",
        help="input .csv file path",
        nargs=1,
        type=str
    )
    args = args_parser.parse_args()
    
    # Read file path
    file_path = args.input[0]

    # Constructure a network using csv file
    graph = {}
    with open(file_path, "r") as f:
        # Read all lines of the file
        lines = f.read().splitlines() 
        for index in tqdm(range(len(lines)), desc="Reading file"):
            line = lines[index]
            node_names = line.split(",")         # split a line into start and stop
            if (len(node_names) != 2):
                # Avoid error
                continue
            # Add nodes to adjacency list
            node1 = int(node_names[0])
            node2 = int(node_names[1])
            # It's a undirected graph, so we add twice
            if (graph.get(node1) is None):
                # If the node is not already in the graph, initialize it
                graph[node1] = np.array([node2])
            else:
                # If already exist, add to it.
                graph[node1] = np.append(graph[node1], node2)
            if (graph.get(node2) is None):
                graph[node2] = np.array([node1])
            else:
                graph[node2] = np.append(graph[node2], node1)
    
    # Number of nodes
    print("Number of nodes: ", len(graph))

    # Number of edges
    degree = np.array([])
    for key in graph:
        degree = np.append(degree, len(graph[key]))
    num_edge = np.sum(degree)
    print("Number of edges: ", int(num_edge/2))

    # Fit power law distribution
    result = powerlaw.Fit(degree)
    result.power_law.xmin

    plt.figure(figsize=(6,5), dpi=400)
    hist, bins = np.histogram(degree, bins=20)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(degree, bins=logbins, density=True, label="Histogramm of degrees")
    x = np.arange(1, np.max(degree),1)
    y = x ** (-result.power_law.alpha)
    plt.plot(x,y,label="Approximates power law function by powerlaw package")
    plt.plot(x, x**(-2), label="Power law with exponent 2")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Degree of node")
    plt.ylabel("Number of nodes")
    plt.legend()
    plt.show()

    ### Find connected components
    sys.setrecursionlimit(len(graph))  # Avoid Recursion Error
    nodes = list(graph)                # Get all node names
    visited = {}                       # Use dictionary get function to speed up, in function is to slow
    cc = []                            # List to save connected components
    for i in tqdm(range(len(nodes)), desc="Searching for connected components"):
        if (visited.get(nodes[i]) is None):
            temp = np.array([])
            temp, visited = DFS(graph, nodes[i], temp, visited)
            cc.append(temp)
    print("Number of connected components: ", len(cc))
    
    # Find the largest connected components
    length = []
    for i in cc:
        length.append(len(i))
    cc_nodes = cc[np.argmax(length)]
    print("The largest connected components contains ", len(cc_nodes), " nodes")

    # Using Dijkstra Algorithm to calculate the one point to all distance
    cc_num = len(cc_nodes)
    num_samples = 100  
    srcs = np.random.choice(np.arange(cc_num), num_samples, replace = False)  # Random choice 100 nodes to compute the distance
    distances = np.array([])
    for i in tqdm(range(num_samples), desc="Calculating distance"):
        # Dijkstra Algorithm
        src = srcs[i]
        distance = np.ones(cc_num) * np.inf
        distance[src] = 0
        queue = [src]
        while (len(queue) != 0):
            current_node = queue.pop(0)
            # Check the neighbors of the current node
            for neighbor in graph[cc_nodes[current_node]]:
                index = np.where(cc_nodes == neighbor)[0][0]
                # Update the distance if it's smaller than the current value
                if (distance[index] > distance[current_node] + 1):
                    distance[index] = distance[current_node] + 1
                    queue.append(index)
        distances = np.append(distances, distance)
    
    # Draw a histogram
    plt.figure(figsize=(12,10), dpi=400)
    plt.hist(distances, density=True)
    plt.xlabel("Number of degrees")
    plt.ylabel("Probability")
    plt.show()



if __name__ == "__main__":
    main()

