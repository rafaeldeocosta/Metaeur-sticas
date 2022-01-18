from igraph import Graph
# from utils import has_node
from igraph import plot

def kruskal(G, E, edge_costs):
    """
        algoritmo que busca uma árvore geradora mínima para um grafo conexo
        com pesos

        https://pt.wikipedia.org/wiki/Algoritmo_de_Kruskal

        Args:
            G - igraph.Graph
            E   - dict - edges
            edge_costs  - dict - cost of edges

    """

    F = [] # Set of Trees

    # creating a forest F where each vertice is a separeted tree
    for v in G.vs:
        t = Graph()
        t.add_vertex(v["name"])
        F.append(t)

    S = E.copy()    # All edges of the graph

    # sorting edge_costs
    edge_costs = dict(sorted(edge_costs.items(), key=lambda item: item[1]))

    # creating and index list of edges sorted by their costs
    edge_index_list = list(edge_costs.keys())

    while S:

        # Removing an edge with minimum cost
        e_index = edge_index_list[0]    # select index of edge with minimum cost
        e = S[e_index]                  # get the edge with e_index: (vi, vj)
        c_e = edge_costs[e_index]       # cost of edge e
        edge_costs.pop(e_index)         # removing e_index  from edge_costs
        edge_index_list.remove(e_index) # removing e_index  from edge_index_list
        S.pop(e_index)                  # removing e_index  from S

        print("Number of tree in F %s" % len(F))
        print("Number of edges in S %s" % len(S))

        # search trees in the forest F where vertices u and v are
        u = e[0]
        v = e[1]

        t_u = None  # tree in F where the vertex u is
        t_v = None  # tree in F where the vertex v is


        # Find out trees that belong to vertex u and v
        for t in F:
            for w in t.vs:  # w is a vertex in the tree t
                if w["name"] == u:
                    t_u = t
                if w["name"] == v:
                    t_v = t

        # discard edges when vertices belong to the same tree
        if t_u == t_v:
            continue

        # Since the edge e  connects two different trees, add it to F,
        # merge them, if they remain a tree

        #
        # check edges of t_u and t_v to create new tree
        #

        # If t_u and t_v dont have edges, create a new tree
        if len(t_u.es) == 0 and len(t_v.es) == 0:
            new_t = Graph()
            new_t.add_vertex(u)
            new_t.add_vertex(v)
            new_t.add_edge(new_t.vs.find(name=u), new_t.vs.find(name=v))

            F.append(new_t)
            F.remove(t_u)
            F.remove(t_v)

        # If only t_v has edges, insert u in t_v
        elif len(t_u.es) == 0 and len(t_v.es) > 0:  # update t_v with u
            t_v.add_vertex(u)
            t_v.add_edge(t_v.vs.find(name=u), t_v.vs.find(name=v))
            F.remove(t_u)

        # If only t_u has edges, insert v in t_u
        elif len(t_u.es) > 0 and len(t_v.es) == 0: # update t_u with v
            t_u.add_vertex(v)
            t_u.add_edge(t_u.vs.find(name=u), t_u.vs.find(name=v))
            F.remove(t_v)

        # Since t_u and t_v have edges, it requires to create a new tree (new_t)
        # with all vertices and edges of t_u and t_v
        else:

            new_t = Graph()
            for w in t_u.vs:
                new_t.add_vertex(w["name"])
            for w in t_v.vs:
                new_t.add_vertex(w["name"])

            # TODO: Temos que pegar o NOME dos vertíces das arestas em T_U e T_V, e fazer find no new_t


            for l in t_u.es:
                source_vertex = t_u.vs[l.source]["name"]
                target_vertex = t_u.vs[l.target]["name"]
                new_t.add_edge(new_t.vs.find(name=source_vertex),
                    new_t.vs.find(name=target_vertex))

            for l in t_v.es:
                source_vertex = t_v.vs[l.source]["name"]
                target_vertex = t_v.vs[l.target]["name"]
                new_t.add_edge(new_t.vs.find(name=source_vertex),
                    new_t.vs.find(name=target_vertex))

            # inserting the edge e that links t_u and t_v
            # print("ligando t_u e t_v com a aresta %s -> %s " % (u,v))
            new_t.add_edge(new_t.vs.find(name=u), new_t.vs.find(name=v))

            # check if new_t remain a tree
            if not new_t.is_tree():
                continue
            else:
                F.append(new_t)
                F.remove(t_u)
                F.remove(t_v)

    T = F[0]
    return T
