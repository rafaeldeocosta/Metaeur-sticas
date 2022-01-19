from utils import create_graph_from, select_sub_graph, calc_pontuacao, vert_premios, tira_grau1
from kruskal import kruskal
from igraph import plot

if __name__ == "__main__":
    print("trabalho de Meta Heuristica: Simulated Annealing")

    arq = "K100.1"

    G, V, vertex_penalties, E, edge_costs = create_graph_from(arq)

    # G.vs["label"] = G.vs["name"]
    # layout = G.layout("lgl")
    # plot(G, layout=layout)

    terminais = vert_premios(G)

    G_terminais, E_terminais, edge_costs_terminais = select_sub_graph(G, terminais)

    # G_terminais.vs["label"] = G_terminais.vs["name"]
    # plot(G_terminais)

    T, T_edges_list = kruskal(G, E, edge_costs)

    # plot graph
    # T.vs["label"] = T.vs["name"]
    # layout = T.layout("lgl")
    # plot(T, layout=layout)

    T_terminais, T_terminais_edges_list = kruskal(G_terminais, E_terminais, edge_costs_terminais)
    print(T_terminais.is_tree())

    T_terminais, T_terminais_edges_list = tira_grau1(T_terminais, terminais)

    T_terminais.vs["label"] = T_terminais.vs["name"]
    layout = T_terminais.layout("lgl")
    plot(T_terminais, layout=layout)

    sol_T = calc_pontuacao(G, T)
    sol_T_terminais = calc_pontuacao(G, T_terminais)

    print('A pontuação de T é: {}'.format(sol_T))
    print('A pontuação de T_terminais é: {}'.format(sol_T_terminais))


