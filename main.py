from igraph import plot
from kruskal import kruskal
from SA import SA
from utils import calc_pontuacao
from utils import create_graph_from
# from utils import get_graph_of_terminals
from utils import get_solution_list
from utils import remove_costly_leafs
from utils import select_sub_graph
from utils import tira_grau1
from utils import vert_premios


if __name__ == "__main__":
    print("trabalho de Meta Heuristica: Simulated Annealing")

    # arq = "instances/our-instances/E/e01.stp-B"
    arq = 'C01-A'
    # G: Original Graph from instance in arq (igraph.Graph)
    # V: Vertices name of G (not vertex id) (list - [Vi, Vj, ...])
    # E: Edges of G (dict - {edge_index:(Vi, Vj)}); edge_index here is different
    #                                               than edge_index in G
    # vertex_penalties: Weight of vertices of G (dict - {vertice name: weight}
    # edge_costs: Cost of edges of G (dict - {edge_index: cost of edge}
    #                           edge_index here is equals in E
    G, V, vertex_penalties, E, edge_costs = create_graph_from(arq)

    #
    # Plot G
    #
    # G.vs["label"] = G.vs["name"]
    # layout = G.layout("lgl")
    # plot(G, layout=layout)

    # ##################################################
    # Generate initial solution with kruskal algorithm #
    ####################################################
    #
    # T: Tree obtained from G, which could be used as Initial Solution of SA
    #       (igraph.Graph)
    # T_edges_list:
    # T, T_edges_list = kruskal(G, E, edge_costs)
    #
    # Plot T
    #
    # T.vs["label"] = T.vs["name"]
    # layout = T.layout("lgl")
    # plot(T, layout=layout)
    #
    # sol_T = calc_pontuacao(G, T)
    # print("Score of T %s" % sol_T)

    # Pruned_T = remove_costly_leafs(T)
    #
    # Plot Pruned_T
    #
    # Pruned_T.vs["label"] = Pruned_T.vs["name"]
    # layout = Pruned_T.layout("lgl")
    # plot(Pruned_T, layout=layout)
    # sol_Pruned_T = calc_pontuacao(G, Pruned_T)
    # print("Score of Pruned_T %s" % sol_Pruned_T)

    # ###################################################################
    # Generate initial solution with subgraph of G connecting terminals #
    #####################################################################

    terminais = vert_premios(G)
    G_terminais, E_terminais, edge_costs_terminais = \
        select_sub_graph(G, terminais)
    #
    # Plot G_Terminais
    #
    # G_terminais.vs["label"] = G_terminais.vs["name"]
    # layout = G.layout("lgl")
    # plot(G_terminais)

    T_terminais, T_terminais_edges_list = \
        kruskal(G_terminais, E_terminais, edge_costs_terminais)

    sol_T_terminais = calc_pontuacao(G, T_terminais)
    print("Score of T_Terminais %s" % sol_T_terminais)

    #
    # Plot T_Terminais
    #
    # T_terminais.vs["label"] = T_terminais.vs["name"]
    # layout = T_terminais.layout("lgl")
    # plot(T_terminais, layout=layout)

    Pruned_T_Terminais = remove_costly_leafs(T_terminais)

    sol_Pruned_T_terminais = calc_pontuacao(G, Pruned_T_Terminais)
    print("Score of Pruned_T_Terminais %s" % sol_Pruned_T_terminais)

    # Pruned_T_Terminais.vs["label"] = Pruned_T_Terminais.vs["name"]
    # plot(Pruned_T_Terminais)


    #
    # Run Simulated Annealing
    #

    #TODO: get params from file called 'cooling_strategies'
    S = Pruned_T_Terminais.copy()
    Temp_ini = 1000
    Temp_fin = 1
    ALPHA = 0.5
    SA_MAX = 400
    f = 'Temp_curr*0.9'  # Função decaimento da temperatura > IMPORTANTE: em função de "Temp_curr"

    Star_S = SA(G, Temp_ini, Temp_fin, S, SA_MAX, f)
    sol_Star_S = calc_pontuacao(G, Star_S)
    Star_S.vs["label"] = Star_S.vs["name"]
    # plot(Star_S)
    print("Score of Star_S %s" % sol_Star_S)
