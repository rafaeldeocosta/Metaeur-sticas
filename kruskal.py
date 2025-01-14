from igraph import Graph
from igraph import plot
from utils import get_atributes

def kruskal(G, E, edge_costs):
    """
        algoritmo que busca uma árvore geradora mínima para um grafo conexo
        com pesos

        https://pt.wikipedia.org/wiki/Algoritmo_de_Kruskal

        Args:
            G - igraph.Graph
            E   - dict - edges
            edge_costs  - dict - cost of edges

        Return:
            T - igraph.Graph - árvore geradora mínima
            edges_list -

    """
    G_ori = G.copy()  # G without alterations

    F = []  # Set of Trees

    # creating a forest F where each vertice is a separeted tree
    for v in G.vs:
        t = Graph()
        t.add_vertex(v["name"])
        F.append(t)

    S = E.copy()    # All edges of the graph

    # sorting edge_costs
    edge_costs = dict(sorted(edge_costs.items(), key=lambda item: item[1]))

    # creating a list of edge_index sorted by their costs
    edge_index_list = list(edge_costs.keys())

    while S:
        #
        # Step #1: Removing an edge with minimum cost
        #
        e_index = edge_index_list[0]    # select index of edge with minimum cost
        e = S[e_index]                  # get the edge with e_index: (vi, vj)
        c_e = edge_costs[e_index]       # cost of edge e
        edge_costs.pop(e_index)         # removing e_index  from edge_costs
        edge_index_list.remove(e_index) # removing e_index  from edge_index_list
        S.pop(e_index)                  # removing e_index  from S


        #
        # Step #2: Selecting vertices u and v of edge e
        #
        u = e[0]
        v = e[1]

        #
        # Step #3: search trees in the forest F that contains vertices u and v
        #

        t_u = None  # tree in F that contain u (igraph.Graph)
        t_v = None  # tree in F that contain v (igraph.Graph)


        # Check trees t in F that contains vertices u and v
        for t in F:
            for w in t.vs:  # w is a vertex in the tree t
                if w["name"] == u:
                    t_u = t
                if w["name"] == v:
                    t_v = t

        # discard edges when vertices belong to the same tree
        if t_u == t_v:
            continue

        #
        #   Step #4: Connect t_u and t_v and update forest F
        #
        #   Since edge e has vertices of different trees, try to merge them,
        #   but if the merged tree stop being a tree, discard such edge
        #

        #
        #   check edges of t_u and t_v to create new tree
        #

        # If t_u and t_v dont have edges, create a new tree with u and v
        # and an edge between them
        if len(t_u.es) == 0 and len(t_v.es) == 0:
            new_t = Graph()
            new_t.add_vertex(u)
            new_t.add_vertex(v)
            new_t.add_edge(new_t.vs.find(name=u), new_t.vs.find(name=v))

            F.append(new_t)
            F.remove(t_u)
            F.remove(t_v)

        # If only t_v has edges, updated t_v with vertice u and create edge e
        # in t_v
        elif len(t_u.es) == 0 and len(t_v.es) > 0:  # update t_v with u
            t_v.add_vertex(u)
            t_v.add_edge(t_v.vs.find(name=u), t_v.vs.find(name=v))
            F.remove(t_u)

        # If only t_u has edges, updated t_u with vertice v and create edge e
        # in t_u
        elif len(t_u.es) > 0 and len(t_v.es) == 0: # update t_u with v
            t_u.add_vertex(v)
            t_u.add_edge(t_u.vs.find(name=u), t_u.vs.find(name=v))
            F.remove(t_v)

        # Since t_u and t_v have edges, it requires to create a new tree (new_t)
        # with all vertices and edges of t_u and t_v. Thus create edge e
        else:

            new_t = Graph()
            for w in t_u.vs:
                new_t.add_vertex(w["name"])
            for w in t_v.vs:
                new_t.add_vertex(w["name"])

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
            new_t.add_edge(new_t.vs.find(name=u), new_t.vs.find(name=v))

            # check if new_t not remain a tree; discard new_t and edge e
            if not new_t.is_tree():
                continue
            else:
                F.append(new_t)
                F.remove(t_u)
                F.remove(t_v)

    # In the end, it is expected that Forest F contains only one tree 
    T = F[0]

    T = get_atributes(G_ori, T)  # Atribui os atributos de vertice e arestas a T

    V_t = T.get_vertex_dataframe().copy()  # dataframe com vertices de aux_G

    V_guide = V_t['name'].to_dict()  # dicionário com {index do vertice em aux_G: Nome do vertice em aux_G}

    T_edges = T.get_edge_dataframe().copy()  # dataframe com arestas de aux_G
    T_edges['source'].replace(V_guide, inplace=True)  # Substitui o valor de source (que e index) para NOME
    T_edges['target'].replace(V_guide, inplace=True)  # Substitui o valor de target (que e index) para NOME

    T_edges['edge_name'] = T_edges[['source', 'target']].apply(tuple, axis=1)  # Cria coluna com as tuplas (source, target)

    edges_list = T_edges['edge_name'].to_list()

    return T, edges_list
