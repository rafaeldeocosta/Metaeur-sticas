from copy   import deepcopy
from random import randint
from igraph import plot
from utils import calc_pontuacao
from random import uniform
from math import exp

def get_neighbor(G, S):
    """
        função que retorna uma solução vizinha de forma aleatória na vizinhança
            de S

            Args:
                G - igraph.Graph - original Graph
                S - igraph.Graph - Tree
    """

    #
    # Plot S
    #
    # S.vs["label"] = S.vs["name"]
    # layout = S.layout("lgl")
    # plot(S, layout=layout)

    S_ = S.copy()

    #
    # First Step: Ramdonly choose number of edges to modify in S
    #

    # number_of_edges_to_change = randint(1, len(S_.es))  #TODO: mudar
    number_of_edges_to_change = 1
    # print("number_of_edges_to_change %s " % number_of_edges_to_change)

    # print(S_.get_vertex_dataframe())
    # print(S_.get_edge_dataframe())
    # print(S_)

    G_orig = G.copy()
    for n_e in range(0,number_of_edges_to_change) :
        G = G_orig.copy()
        #
        #   Step #2: Ramdonly choose an edge to remove
        #
        e_to_remove = S_.es[randint(0, len(S_.es) - 1)]
        # print(e_to_remove)
        S_.delete_edges(e_to_remove)
        G.delete_edges(e_to_remove)

        #
        # Step #3: Select trees generated after removing such edge
        #
        clusters = S_.clusters()
        if len(clusters) == 2:
            t_u = clusters.subgraph(0)
            t_v = clusters.subgraph(1)
        else:
            print("there is no two trees after removing an edge")

        #
        # Step #4: ramdonly choose vertice u from t_u and v from t_v
        #
        u = t_u.vs[randint(0, len(t_u.vs) - 1)]["name"]
        v = t_v.vs[randint(0, len(t_v.vs) - 1)]["name"]

        # print("u = %s" % u)
        # print("v = %s" % v)

        #
        # Step #5: Get in the original Graph (G) all paths that connects u and v
        #
        P = G.get_all_shortest_paths(G.vs.find(name=u),
                                        to=G.vs.find(name=v))
        if len(P) == 0:
            S_ = S.copy()
            continue

        #
        # Step #6: ramdonly choose one path from P
        #
        p = P[randint(0, len(P)-1)] # list of vertex_id of G
        # print("caminho escolhido")
        # print(p)
        # for x in p:
        #     print(G.vs[x])

        #
        # Step #7: Trying to reconnect the tree S_ using p
        #

        # tenta reconectar S_ (as arvores t_u e t_v)
        # se não gerar ciclo, aceita como vizinho

        # print("=== S_ before ===")
        # print(S_)
        #
        # t_u.vs["label"] = t_u.vs["name"]
        # layout = t_u.layout("lgl")
        # plot(t_u, layout=layout)
        #
        # t_v.vs["label"] = t_v.vs["name"]
        # layout = t_v.layout("lgl")
        # plot(t_v, layout=layout)

        s = deepcopy(u)
        for t in p:
            vertex_t_in_G = G.vs[t]["name"]

            if  s == vertex_t_in_G:
                # discarding the first element of p because it is u
                continue

            else:

                #
                # Step #8: Checking if vertices s and t are in S_
                #               since t is a vertice name in G and may not
                #               exist in S_ yet

                # getting vertice name of vertex id t in G
                if len(S_.vs.select(name=s)) == 0:
                    S_.add_vertex(name=s)
                    vertex_s = S_.vs.find(name=s)

                    if vertex_s['penalties'] is None:
                        vertex_s['penalties'] = G.vs.find(name=s)['penalties']

                elif len(S_.vs.select(name=vertex_t_in_G)) == 0:
                    S_.add_vertex(name=vertex_t_in_G)
                    vertex_t = S_.vs.find(name=vertex_t_in_G)

                    if vertex_t['penalties'] is None:
                        vertex_t['penalties'] = \
                            G.vs.find(name=vertex_t_in_G)['penalties']

                #
                # Step #9: Before add the edge s --t , check if the edge
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

                else:
                    s = deepcopy(vertex_t_in_G)
                    continue


            s = deepcopy(vertex_t_in_G)

        #
        # Step 10: S_ after reconnection should be a tree.
        #
        if not S_.is_tree():
            S_ = S.copy()

        # print("=== S_ after ===")
        # print(S_)
        # print(S_.get_vertex_dataframe())
        # print(S_.get_edge_dataframe())

        # S_.vs["label"] = S_.vs["name"]
        # layout = S_.layout("lgl")
        # plot(S_, layout=layout)

    return S_


def metropolis(G, S_ini, Temperature, n_iter_temperature):
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
        S_viz = get_neighbor(G, Best_S)

        if  calc_pontuacao(G, S_viz) < calc_pontuacao(G, Best_S):
            Best_S = S_viz.copy()
        else:
            DELTA =  calc_pontuacao(G, S_viz) - calc_pontuacao(G, Best_S)
            p = uniform(0,1)
            # print("p = %s" % p)
            # print("Temperature = %s" % Temperature)
            # print("DELTA = %s" % DELTA)
            # print("e = %s" % exp(DELTA/Temperature))
            if  exp(DELTA/Temperature) > 0:
                if p < (1 / exp(DELTA/Temperature)):
                    Best_S = S_viz.copy()

    return Best_S
