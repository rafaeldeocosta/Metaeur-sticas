from copy import deepcopy
from random import randint
from random import choice
from igraph import plot
from utils import calc_pontuacao
from random import uniform
from math import exp
from math import inf

import igraph as ig

def get_neighbor(G, S, T_curr, Temp_ini,points=[1, 1, 1]):

    # Step  # 1: Choose which method to use
    methods = ['remove_vertices' for item in range(0, points[0])] + \
              ['remove_edges' for item in range(0, points[1])] + \
              ['both' for item in range(0, points[2])]

    process = choice(methods)

    # Step  # 2: Choose the number of changes
    # TODO: NUMERO DE MOFICAÇÕES PROPORCIONAL A TEMPERATURA
    max_changes = round(len(S.vs)*(T_curr/Temp_ini))
    max_changes = max_changes if max_changes > 1 else 1

    n_changes = choice(list(range(1, max_changes+1)))

    # Step  # 2: Get neighbor with the choosen process
    if process == 'remove_vertices':
        if n_changes < len(S.vs) - 2:
            S_ = remove_vertices(G, S, n_changes)
        else:
            S_ = S.copy()

    if process == 'remove_edges':
        S_ = remove_edges(G, S, n_changes)

    if process == 'both':

        # TODO: ESCOLHER NUMERO DE MODIFICADOES PARA CADA UM
        S_ = S.copy()

        for n in range(0, n_changes):
            next = choice(['remove_vertices', 'remove_edges'])

            if next == 'remove_vertices':

                if n_changes < len(S_.vs) - 2:
                    S_ = remove_vertices(G, S_, 1)
                else:
                    continue

            if next == 'remove_edges':
                S_ = remove_edges(G, S_, 1)

    return S_, process


def join_forest(G, S_):

    """
        Function that returns a single tree made from a forest.

    :param G: igraph.Graph
        The original Graph.
    :param S_: igraph.Graph
        A forest.
    :return: igraph.Graph
        A tree.
    """

    # Step  # 1: Genarate list of trees in S_ forest
    forest = S_.clusters().subgraphs()  # List of trees in S_

    while len(forest) > 1:

        # print(forest)
        # Step  #2: Select trees generated after removing such edge
        t_u = forest[0]
        t_v = forest[1]

        # Step #3: ramdonly choose vertice u from t_u and v from t_v
        u = choice(t_u.vs)['name']
        v = choice(t_v.vs)['name']

        # Step #4: Get in the original Graph (G) all paths that connects u and v
        paths = G.get_all_shortest_paths(G.vs.find(name=u),
                                         to=G.vs.find(name=v))

        # Step #5: ramdonly choose one path from P TODO: THREADS
        if paths:
            p = choice(paths)
        else:
            print('entrou no continue do paths')
            continue

        # First s is u
        s = u
        for t in p:

            vertex_t_in_G = G.vs[t]["name"]

            if s == vertex_t_in_G:
                # discarding the first element of p because it is u
                continue

            else:

                # Step #6: Checking if vertices s and t are in S_
                #               since t is a vertice name in G and may not
                #               exist in S_ yet

                # getting vertice name of vertex id t in G
                if s not in S_.vs['name']:  # len(S_.vs.select(name=s)) == 0:
                    S_.add_vertex(name=s)
                    vertex_s = S_.vs.find(name=s)

                    if vertex_s['penalties'] is None:
                        vertex_s['penalties'] = G.vs.find(name=s)['penalties']

                elif vertex_t_in_G not in S_.vs['name']:  # len(S_.vs.select(name=vertex_t_in_G)) == 0:
                    S_.add_vertex(name=vertex_t_in_G)
                    vertex_t = S_.vs.find(name=vertex_t_in_G)

                    if vertex_t['penalties'] is None:
                        vertex_t['penalties'] = \
                            G.vs.find(name=vertex_t_in_G)['penalties']

                #
                # Step #7: Before add the edge s --t , check if the edge
                #           between s and t already exist
                #
                if len(S_.es.select(_source=S_.vs.find(name=s),
                                    _target=S_.vs.find(name=vertex_t_in_G))) == 0 and \
                        len(S_.es.select(_source=S_.vs.find(name=vertex_t_in_G),
                                         _target=S_.vs.find(name=s))) == 0:

                    S_.add_edges([(S_.vs.find(name=s),
                                   S_.vs.find(name=vertex_t_in_G))])

                    e_in_S_ = S_.es.select(_source=S_.vs.find(name=s),
                                           _target=S_.vs.find(name=vertex_t_in_G))

                    e_in_G = G.es.select(_source=G.vs.find(name=s),
                                         _target=G.vs.find(name=vertex_t_in_G))

                    # Updating cost of new edge in S_
                    if len(e_in_G) == 0:
                        print("edge does not exist in G")
                    else:
                        e_in_S_['cost'] = e_in_G[0]['cost']

                    # TODO: MARCAR PARA MOSTRAR PRO RAFAEL
                    # plot = S_.copy()
                    # plot.vs["label"] = plot.vs["name"]
                    # ig.plot(plot).show()

                    #
                    # Step #7: To verify in which tree edges are been added, to check
                    #               if has an cycle
                    #
                    for subgraph in S_.clusters().subgraphs():
                        if s in subgraph.vs['name']:
                            tree = subgraph.copy()
                            break
                        else:
                            print('Join_forest > s not in subgraph')
                            tree = []

                    # if when e_in_S_ was added a cicle was created, dele edege
                    if tree and not tree.is_tree():
                        S_.delete_edges(e_in_S_)

                        # TODO: MARCAR PARA MOSTRAR PRO RAFAEL
                        # plot = S_.copy()
                        # plot.vs["label"] = plot.vs["name"]
                        # ig.plot(plot).show()

                    # If S_ is tree, can stop loop < Stopping criterion
                    if S_.is_tree():
                        break

                else:
                    s = vertex_t_in_G
                    continue

            s = vertex_t_in_G

        # List of forests
        forest = S_.clusters().subgraphs()

    return S_

def remove_edges(G, S, n_changes=1):
    """
        Dado um número de alterações, remove arestas e reconecta as árvores formadas,
        retornando no final uma solução vizinha S_
    :param G_aux: igraph.Graph
        Original graph
    :param S: igraph.Graph
        A tree of solution
    :param n_changes: int
        Number of changes.
    :return: igraph.Graph
        New solution.
    """

    S_ = S.copy()
    G_aux = G.copy()
    G_orig = G_aux.copy()

    for n_e in range(0, n_changes):

        G_aux = G_orig.copy()  # TODO: Acho que não tem que ter isso

        # Step #1: Ramdonly choose an edge to remove
        e_to_remove = choice(S_.es)

        source = e_to_remove.tuple[0]  # edge to remove source index in S_
        target = e_to_remove.tuple[1]  # edge to remove target index in S_

        source = S_.vs.find(source)['name']  # edge to remove source NAME in S_
        target = S_.vs.find(target)['name']  # edge to remove target NAME in S_

        # Edge to remove index in G
        e_to_remove_in_G = G_aux.es.select(_source=G_aux.vs.find(name=source),
                                           _target=G_aux.vs.find(name=target)).indices[0]

        G_aux.delete_edges(e_to_remove_in_G)

        S_.delete_edges(e_to_remove)

        # If S_ is tree, keep up with alterations < will never enter here
        if S_.is_tree():
            continue

        # if S_ is a forest, join trees in forest
        else:
            S_ = join_forest(G_aux, S_)

    return S_

def remove_vertices(G, S, n_changes=1):
    """
        Dado um número de alterações, remove vértices e reconecta as árvores formadas,
        retornando no final uma solução vizinha S_
    :param G_aux: igraph.Graph
        Original graph
    :param S: igraph.Graph
        A tree of solution
    :param n_changes: int
        Number of changes.
    :return: igraph.Graph
        New solution.
    """
    S_ = S.copy()
    G_aux = G.copy()
    G_orig = G_aux.copy()

    for n in range(0, n_changes):

        G_aux = G_orig.copy()

        # print(S_)
        v_to_remove = choice(S_.vs)  # random vertice in S_

        edges_from_v_S = S_.es.select(_source=v_to_remove).indices  # indexes of edges from v_to_remove in S_

        edges_to_remove_G = []

        # Get index of edges to remove in G
        for edge in edges_from_v_S:

            idx_u = S_.es.find(edge).tuple[0]  # index of source from edge in S_
            idx_v = S_.es.find(edge).tuple[1]  # index of target from edge in S_

            u = S_.vs.find(idx_u)['name']  # name of source from edge
            v = S_.vs.find(idx_v)['name']  # name of target from edge

            edge_in_G = G_aux.es.select(_source=G_aux.vs.find(name=u), _target=G_aux.vs.find(name=v)).indices[0]

            edges_to_remove_G.append(edge_in_G)

        G_aux.delete_edges(edges_to_remove_G)

        #
        v_to_remove_G = G_aux.vs.find(name=v_to_remove['name'])  # index of v_to_remove in G

        G_aux.delete_vertices(v_to_remove_G)

        S_.delete_vertices(v_to_remove)

        # If S_ is tree, keep up with alterations
        if S_.is_tree():
            continue

        # if S_ is a forest, join trees in forest
        else:
            S_ = join_forest(G_aux, S_)

    return S_

def metropolis(G, S_ini, Temperature, Temp_ini,n_iter_temperature):
    """
        função que implementa o algoritmo de metropolis

        Args
            G - igraph.Graph - original Graph
            S_ini - igraph.Graph - initial Solution of PCSTP ("Steiner" Tree)
            Temperature - float - Temperature used to calculate probability to
                accept worst solutions
            n_iter_temperature - maximum of iteration in this Temperature

        Returns
            S_ - igraph.Graph - Solution of PCSTP ("Steiner" Tree)

    """

    print("Running metropolis for temperature %s" % Temperature)

    Best_S = S_ini.copy()

    for i in range(0, n_iter_temperature):

        # GERAR um vizinho s de forma aleatória na vizinhança ℵ𝑠
        S_viz, process = get_neighbor(G, Best_S, Temperature, Temp_ini)

        # S_viz = remove_vertices(G, Best_S)

        if calc_pontuacao(G, S_viz) <= calc_pontuacao(G, Best_S):
            Best_S = S_viz.copy()
        else:
            DELTA = calc_pontuacao(G, S_viz) - calc_pontuacao(G, Best_S)
            p = uniform(0, 1)
            # print("p = %s" % p)
            # print("Temperature = %s" % Temperature)
            # print("DELTA = %s" % DELTA)
            # print("e = %s" % exp(DELTA/Temperature))

            try:
                ans = exp(DELTA/Temperature)
            except OverflowError:
                ans = inf

            if ans > 0 and ans != inf:
                if p < (1 / ans):
                    Best_S = S_viz.copy()

    return Best_S
