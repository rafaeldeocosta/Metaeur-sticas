from utils import create_graph_from
from kruskal import kruskal
from igraph import plot

if __name__ == "__main__":
    print("trabalho de Meta Heuristica: Simulated Annealing")

    arq = "K100.1"

    G, V, vertex_penalties, E, edge_costs = create_graph_from(arq)
    terminais = [4, 6, 13, 17, 43, 46, 61, 71, 78, 79, 91]

    T = kruskal(G, E, edge_costs)

    # plot graph
    T.vs["label"] = T.vs["name"]
    layout = T.layout("lgl")
    plot(T, layout=layout)


