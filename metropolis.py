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

    S_ = S.copy()

    # S_.vs["label"] = S_.vs["name"]
    # layout = S_.layout("lgl")
    # plot(S_, layout=layout)


    # escolhe aleatoriamente a quantidade de arestas modificadas da arvore
    #   pra remover
    number_of_edges_to_change = randint(1, len(S_.es))  #TODO: mudar
    number_of_edges_to_change = 1
    # print("number_of_edges_to_change %s " % number_of_edges_to_change)

    # print(S_.get_vertex_dataframe())
    # print(S_.get_edge_dataframe())
    # print(S_)

    n_e = 1
    while n_e <= number_of_edges_to_change:
        # escolhe aleatoriamente vertices u e v das arvores geradas
        #   pela remoção da aresta

        e_to_remove = S_.es[randint(0, len(S_.es) - 1)]
        # print(e_to_remove)
        S_.delete_edges(e_to_remove)

        # selecting generated trees after removing an edge

        clusters = S_.clusters()
        if len(clusters) == 2:
            t_u = clusters.subgraph(0)
            t_v = clusters.subgraph(1)
        else:
            print("there is no two tree after removing an edge")

        # escolhe os vertices u e v das duas arvores criadas após a
        #   remoção da aresta acima

        u = t_u.vs[randint(0, len(t_u.vs) - 1)]["name"]
        v = t_v.vs[randint(0, len(t_v.vs) - 1)]["name"]

        # print("u = %s" % u)
        # print("v = %s" % v)

        # roda dijkstra entre u e v (se não me engano, o dijkstra
        #   retorna varios caminhos)

        P = G.get_all_shortest_paths(G.vs.find(name=u), to=G.vs.find(name=v))

        if len(P) == 0:
            continue

        # escolhe aleatorio dentre os caminhos de dijkstra

        p = P[randint(0, len(P)-1)] # list of vertex_id of G
        # print("caminho escolhido")
        # print(p)
        # for x in p:
        #     print(G.vs[x])

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

        s = u
        for t in p:
            if G.vs[t]["name"] == s:
                continue

            else:
                if len(S_.vs.select(name=s)) == 0:
                    S_.add_vertex(name=s)
                    s_s = S_.vs.find(name=s)
                    if s_s['penalties'] is None:
                        s_s['penalties'] = G.vs.find(name=s)['penalties']


                elif len(S_.vs.select(name=G.vs[t]["name"])) == 0:
                    S_.add_vertex(name=G.vs[t]["name"])
                    s_t = S_.vs.find(name=G.vs[t]["name"])

                    if s_t['penalties'] is None:
                        s_t['penalties'] = G.vs.find(name=G.vs[t]["name"])['penalties']



                # Before add edge, check if the edge between s and t exist
                if len(S_.es.select(_source=S_.vs.find(name=s),
                                    _target=S_.vs.find(name=G.vs[t]["name"]))) == 0 and \
                   len(S_.es.select(_source=S_.vs.find(name=G.vs[t]["name"]),
                                    _target=S_.vs.find(name=s))) == 0:

                    S_.add_edges([(S_.vs.find(name=s),
                                    S_.vs.find(name=G.vs[t]["name"]))])


                    e = S_.es.select(_source=S_.vs.find(name=s),
                                    _target=S_.vs.find(name=G.vs[t]["name"]))
                    e_g = G.es.select(_source=G.vs.find(name=s),
                                    _target=G.vs.find(name=G.vs[t]["name"]))

                    e['cost'] = e_g[0]['cost']

                s = G.vs[t]["name"]

        # print("=== S_ after ===")
        # print(S_)
        # print(S_.get_vertex_dataframe())
        # print(S_.get_edge_dataframe())

        # S_.vs["label"] = S_.vs["name"]
        # layout = S_.layout("lgl")
        # plot(S_, layout=layout)

        # senão, volta e escolhe outro caminho de dijkstra

        n_e += 1

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

    print("Running metropolis for Temperature %s" % Temperature)

    Best_S = S_ini

    i = 1
    while i < n_iter_temperature:

        # GERAR um vizinho s de forma aleatória na vizinhança ℵ𝑠
        S_viz = get_neighbor(G, Best_S)

        if calc_pontuacao(G, Best_S) < calc_pontuacao(G, S_viz):
            Best_S = S_viz
        else:
            DELTA =  calc_pontuacao(G, S_viz) - calc_pontuacao(G, Best_S)
            p = uniform(0,1)
            if p < (1 / exp(DELTA/Temperature)):
                Best_S = S_viz
        i += 1

    return Best_S
