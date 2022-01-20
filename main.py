from igraph import plot
from kruskal import kruskal
from utils import calc_pontuacao
from utils import create_graph_from
from utils import remove_costly_leafs
from utils import select_sub_graph
from utils import tira_grau1
from utils import vert_premios


if __name__ == "__main__":
    print("trabalho de Meta Heuristica: Simulated Annealing")

    arq = "K100.1"

    G, V, vertex_penalties, E, edge_costs = create_graph_from(arq)

    # G.vs["label"] = G.vs["name"]
    # layout = G.layout("lgl")
    # plot(G, layout=layout)

    # G_terminais.vs["label"] = G_terminais.vs["name"]
    # plot(G_terminais)

    T, T_edges_list = kruskal(G, E, edge_costs)

    # plot graph
    # T.vs["label"] = T.vs["name"]
    # layout = T.layout("lgl")
    # plot(T, layout=layout)

    terminais = vert_premios(G)
    G_terminais, E_terminais, edge_costs_terminais = select_sub_graph(G, terminais)

    T_terminais, T_terminais_edges_list = kruskal(G_terminais, E_terminais, edge_costs_terminais)
    # print(T_terminais.is_tree())

    T_terminais, T_terminais_edges_list = tira_grau1(T_terminais, terminais)

    T_terminais.vs["label"] = T_terminais.vs["name"]
    layout = T_terminais.layout("lgl")
    plot(T_terminais, layout=layout)


    Pruned_T = remove_costly_leafs(T_terminais)

    # Pruned_T.vs["label"] = Pruned_T.vs["name"]
    # layout = Pruned_T.layout("lgl")
    # plot(Pruned_T, layout=layout)

    sol_T = calc_pontuacao(G, T)
    sol_T_terminais = calc_pontuacao(G, T_terminais)
    sol_Pruned_T = calc_pontuacao(G, Pruned_T)

    print('A pontuação de T é: {}'.format(sol_T))
    print('A pontuação de T_terminais é: {}'.format(sol_T_terminais))
    print('A pontuação de Pruned_T é: {}'.format(sol_Pruned_T))
