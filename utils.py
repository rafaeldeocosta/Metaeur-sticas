from igraph import Graph
from igraph import plot

def create_graph_from(f):
    """
        create graph from file whose format is like K100.1 file

        Args:
            f - string - filename of the PCSTP instance

        Returns:
            G: igraph.Graph
            V:  list  of vertices
            vertex_penalties: dict - {vertice: weight of vertices}
            E: dict - {index: edges (Vi, Vj)}
            edge_costs: dict - {index: cost of edges}

    """

    fd = open(f, "r")
    instance = fd.readlines()

    node_index = instance.index("node\n")
    link_index = instance.index("link\n")


    V = []  # list of vertices [v1,v2,...)]
    vertex_penalties = {}   # dict of penalties which key is the vertex id
                            # and value is the weight of that vertice
                            # {v1:w1, v2:w2, ...}

    for v in instance[node_index+2:link_index]:
        v = v.split()
        V.append(int(v[0]))   # vertices ids are subtracted by one because
                                # igraph vertex id begins with 0
        vertex_penalties[int(v[0])] = int(v[3])


    E = {}  # dict of edges which key is the edge id and value the edge (vi,vj)
            # {e1:(vi, vj), e2:(vj, vk), ...}

    edge_costs = {} # dict of edge costs which key is the edge id and value is
                    # the cost of that edge {e1:c1, e2:c2, ...}


    for e in instance[link_index+2:-1]:
        e = e.split()
        E[int(e[0])] = ( int(e[1])-1, int(e[2])-1)
        edge_costs[int(e[0])] = int(e[3])


    G = Graph()
    G.add_vertices(len(V))  # adding the number of vertices of the graph
    G.vs["name"] = list(V) # setting labels to identify vertices
    G.add_edges(list(E.values()))  # inserting edges of the graph

    G.es['weight'] = list(edge_costs.values())
    G.vs['cost'] = list(vertex_penalties.values())

    return G, V, vertex_penalties, E, edge_costs

if __name__ == "__main__":
    arq = "K100.1"
    G, V, vertex_penalties, E, edge_costs = create_graph_from(arq)
    print(V)
    print(vertex_penalties)
    print(E)
    print(edge_costs)

    layout = G.layout("lgl")
    plot(G,layout=layout)
