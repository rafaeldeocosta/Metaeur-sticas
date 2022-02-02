import os

import pandas as pd

from igraph import Graph
from igraph import plot
import itertools as it



def create_graph_from_stp(f):
    """
        create graph from file whose format is stp

        Args:
            f - string - filename of the PCSTP instance in stp format

        Returns:
            G: igraph.Graph
            V:  list  of vertices
            vertex_penalties: dict - {vertice: weight of vertices}
            E: dict - {index: edges (Vi, Vj)}
            edge_costs: dict - {index: cost of edges}

    """

    fd = open(f, "r")
    instance = fd.readlines()
    # for i in instance:
    #     print(i)
    # exit()

    graph_section_index = instance.index("SECTION Graph\n")
    nodes = instance[graph_section_index+1].split()[1]
    V = list(range(1, int(nodes)+1)) # list of vertices [v1,v2,...)]
    vertex_penalties = {}   # dict of penalties which key is the vertex id
                            # and value is the weight of that vertice
                            # {v1:w1, v2:w2, ...}
    for v in V:
        vertex_penalties[v] = 0



    E = {}  # dict of edges which key is the edge id and value the edge (vi,vj)
            # {e1:(vi, vj), e2:(vj, vk), ...}

    edge_costs = {} # dict of edge costs which key is the edge id and value is
                    # the cost of that edge {e1:c1, e2:c2, ...}

    links = instance[graph_section_index+2].split()[1]
    i=1
    for e in instance[graph_section_index+3:graph_section_index+int(links)+3]:
        e = e.split()

        E[i] = ( int(e[1]) , int(e[2]) )
        edge_costs[i] = float(e[3])
        i+=1


    terminals_section_index = instance.index("SECTION Terminals\n")
    terminals = int(instance[terminals_section_index+1].split()[1])
    for t in instance[terminals_section_index+2:
                terminals_section_index+int(terminals)+2]:
        t = t.split()
        vertex_penalties[int(t[1])] = float(t[2])

    G = Graph()
    G.add_vertices(len(V))  # adding the number of vertices of the graph
    G.vs["name"] = list(V) # setting labels to identify vertices

    edges_list = []
    for e in E.values():
        edges_list.append((e[0]-1, e[1]-1))

    # inserting edges of the graph
    G.add_edges(edges_list)

    # inserting edges of the graph
    # for e in E.values():
    #     u = e[0]
    #     v = e[1]
    #     G.add_edge(G.vs.find(name=u), G.vs.find(name=v))

    # Updating G with cost of edges and penalties of vertices
    G.es['cost'] = list(edge_costs.values())
    G.vs['penalties'] = list(vertex_penalties.values())

    return G, V, vertex_penalties, E, edge_costs


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
        vertex_penalties[int(v[0])] = float(v[3])


    E = {}  # dict of edges which key is the edge id and value the edge (vi,vj)
            # {e1:(vi, vj), e2:(vj, vk), ...}

    edge_costs = {} # dict of edge costs which key is the edge id and value is
                    # the cost of that edge {e1:c1, e2:c2, ...}


    for e in instance[link_index+2:-1]:
        e = e.split()
        E[int(e[0])] = ( int(e[1]), int(e[2]))
        edge_costs[int(e[0])] = float(e[3])


    G = Graph()
    G.add_vertices(len(V))  # adding the number of vertices of the graph
    G.vs["name"] = list(V) # setting labels to identify vertices

    # inserting edges of the graph
    for e in E.values():
        u = e[0]
        v = e[1]
        G.add_edge(G.vs.find(name=u), G.vs.find(name=v))

    # Updating G with cost of edges and penalties of vertices
    G.es['cost'] = list(edge_costs.values())
    G.vs['penalties'] = list(vertex_penalties.values())

    return G, V, vertex_penalties, E, edge_costs


def select_sub_graph(G, vertices):

    """
        Função que retorna um subgrafo com os mínimos caminhos entre os
            vértices da lista "vertices"

        :param G: igraph.Graph
            Grafo principal, onde o custo das arestas está no atributo "cost" e
                as penalidades dos véritces em "penalties"
        :param vertices: list
            Lista de NOMES de vertices

        :return:
            aux_G: igraph.Graph
                Subgrafo de G com vértices da lista "vertices"
            E_aux: dict
                Dicionario com {index da arestas: (u, v)} onde u e v são nomes
                    dos vértices source e target
            edges_cost_aux: dict
                Dicionario com {index arestas: cost aresta}
    """

    # combinacao 2 a 2 dos itens na lista de vertices
    combi_V = list(it.combinations(vertices, 2))

    verts_G = G.get_vertex_dataframe().copy()

    aux_G = Graph()
    aux_G.vs['name'] = []

    for item in combi_V:

        source = item[0]
        destination = item[1]

        # retorna os indexes dos vértices do menor caminho
        minimal_path = G.get_shortest_paths(G.vs.find(name=source),
                                            G.vs.find(name=destination),
                                            weights='cost')[0]

        for id, vert in enumerate(minimal_path[:-1]):

            G_u_idx = vert  # Index de u em G
            u = G.vs.find(vert)['name']  # Nome de u

            G_v_idx = minimal_path[id+1]  # Index de v em G
            v = G.vs.find(G_v_idx)['name']  # Nome de v

            if u not in aux_G.vs['name']:

                # penalidade de u
                penalty_u = verts_G.loc[verts_G['name'] == u,
                                                'penalties'].values[0]

                # dessa maneira a penalidade já entra no G
                aux_G.add_vertices([u], attributes={'penalties': [penalty_u]})

            if v not in aux_G.vs['name']:
                # aux_G.add_vertex(v)
                penalty_v = verts_G.loc[verts_G['name'] == v, 'penalties'].values[0]  # penalidade de v
                aux_G.add_vertices([v], attributes={'penalties': [penalty_v]})  # dessa maneira a penalidade já entra no G

            edge_cost = G.es.select(_source=G_u_idx, _target=G_v_idx)['cost'][0]   # custo da aresta u v

            edges_list = aux_G.get_edgelist()   # lista de arestas com INDEX dos vertices
            edges_list = [set(item) for item in edges_list]  # transforma lista de tuplas, para lista de sets

            aux_u_idx = aux_G.vs.find(name=u)  # index de u em aux_G
            aux_v_idx = aux_G.vs.find(name=v)  #index de v em aux_G

            if {aux_u_idx.index, aux_v_idx.index} not in list(edges_list):
                aux_G.add_edges([(aux_u_idx, aux_v_idx)], attributes={'cost': [edge_cost]})  # dessa maneira o custo já entra no G

    V_aux = aux_G.get_vertex_dataframe().copy()  # dataframe com vertices de aux_G

    V_guide = V_aux['name'].to_dict()  # dicionário com {index do vertice em aux_G: Nome do vertice em aux_G}

    edges_aux = aux_G.get_edge_dataframe().copy()  # dataframe com arestas de aux_G
    edges_aux['source'].replace(V_guide, inplace=True)  # Substitui o valor de source (que e index) para NOME
    edges_aux['target'].replace(V_guide, inplace=True)  # Substitui o valor de target (que e index) para NOME

    edges_aux['edge_name'] = edges_aux[['source', 'target']].apply(tuple, axis=1)  # Cria coluna com as tuplas (source, target)

    E_aux = edges_aux['edge_name'].to_dict()  # dicionario com {index da arestas: (u, v)}
    edges_cost_aux = edges_aux['cost'].to_dict()  # dicionario com {index arestas: cost aresta}

    return aux_G, E_aux, edges_cost_aux


def get_atributes(G, tree):
    """
    Atribui os custos de arestas e penalidades de vertices de G em tree

    :param G: igraph.graph
        G e o grafo inicial do problema
    :param tree:
        tree da solucao de kruskal
    :return:
        tree com custos e penalidades
    """

    for v in tree.vs:
        v['penalties'] = G.vs.find(name=v['name'])['penalties'] # procura a penalidade do vertice v em G e insere no vertice v em tree

    for edge in tree.es:

        u = tree.vs.find(edge.tuple[0])['name']
        v = tree.vs.find(edge.tuple[1])['name']

        edge['cost'] = G.es.select(_source=G.vs.find(name=u), _target=G.vs.find(name=v))['cost'][0] # procura a custo da aresta (u,v) em G e insere na aresta (u,v) em tree

    return tree


def calc_pontuacao(G, tree):
    """
    Calcula soma dos custos da arestas com a penalidade dos vértices não escolhidos
    :param G: igraph.graph
        Grafo G original do problema
    :return: float
        Soma dos custos das arestas com a penalidade dos vértices não escolhidos

    """

    vertices_tree = tree.get_vertex_dataframe() # dataframe com vertices de tree

    edges_tree = tree.get_edge_dataframe()  # dataframe de arestas de tree

    vertices_G = G.get_vertex_dataframe()  # dataframe com vertices e penalidades de G

    list_v_tree = vertices_tree['name'].to_list()  # lista de vertices de tree
    list_v_G = vertices_G['name'].to_list()  # lista de vertices de G

    vert_n_escolhido = [item for item in list_v_G if item not in list_v_tree]  # lista de vertices não escolhidos

    soma_penalties = vertices_G.loc[vertices_G['name'].isin(vert_n_escolhido), 'penalties'].sum()  # seleciona a coluna penalidades, do vertices que estao em vert_n_escolhido e soma

    soma_costs = edges_tree['cost'].sum()  # soma da coluna de custos das arestas de tree

    total = soma_penalties + soma_costs

    return total


def tira_grau1(tree, terminais, iter_max=100):
    """
        Looping que retira os vertices de grau 1, até não ter vertices de grau 1 em vertices de stainer.
    :param tree: graph.igraph
        Um grafo para se retirar os vertices de grau1
    :param terminais: list
        Lista de "terminais", vertices que tem peso
    :param iter_max:
        Numero máximo de iterações o padrão é 100

    :return:
        tree: igraph.graph
            tree sem vertices de stainer de grau1
        edges_list: list
            lista de arestas com NOMES (u,v)
    """

    for item in list(range(0, iter_max)):

        vert_grau1 = [index for index, item in enumerate(tree.vs.degree()) if item == 1]  # pega o index dos vertices de grau1

        tree_vertices = tree.get_vertex_dataframe()  # dataframe de vertices

        name_grau1 = tree_vertices.loc[vert_grau1, ['name']]  # seleciona os vertices de grau 1 em tree_vertices

        name_grau1 = name_grau1.loc[~name_grau1['name'].isin(terminais)]  # verifica os vertices de grau 1 que nao sao terminais

        v_retirar = name_grau1.index.to_list()  # pega lista index dos vertices de grau 1 para retirar

        # se tiver vertices para retirar retire se nao para o laço
        if v_retirar:
            tree.delete_vertices(v_retirar)

        else:
            break


    V_tree = tree.get_vertex_dataframe().copy()  # dataframe com vertices de aux_G

    V_guide = V_tree['name'].to_dict()  # dicionário com {index do vertice em aux_G: Nome do vertice em aux_G}

    tree_edges = tree.get_edge_dataframe().copy()  # dataframe com arestas de aux_G
    tree_edges['source'].replace(V_guide, inplace=True)  # Substitui o valor de source (que e index) para NOME
    tree_edges['target'].replace(V_guide, inplace=True)  # Substitui o valor de target (que e index) para NOME

    tree_edges['edge_name'] = tree_edges[['source', 'target']].apply(tuple, axis=1)  # Cria coluna com as tuplas (source, target)

    edges_list = tree_edges['edge_name'].to_list()  # Lista de tuplas com as arestas

    return tree, edges_list


def vert_premios(G):
    """
    Retorna lista de vertices com premios.
    :param G: igraph.Graph
    :return: list

    """
    g_vertices = G.get_vertex_dataframe()
    v_premios = g_vertices.loc[g_vertices['penalties'] != 0, 'name'].to_list()

    return v_premios


def remove_costly_leafs(T):
    """
        função que remove os vértices folha do grafo cujo prize é menor que
            o custo da aresta para conectá-lo

            Args:
                T - igraph.Graph - Tree
                T_terminais_edges_list - list of edges [(i,j),...]

            Return:
                Pruned_T - igraph.Graph - Tree without costly leafs

    """

    Pruned_T = T.copy()

    while True:
        # select leafs, i.e. vertices with degree = 1
        leafs_vertices = []
        for v in Pruned_T.vs:
            if v.degree() == 1:
                leafs_vertices.append(v)

        vertices_to_remove = []
        for v in leafs_vertices:
            # print("vertice %s" % v["name"])
            # print(v)
            for e in v.all_edges():
                # print(e)
                # print(Pruned_T.vs[e.source]["name"])
                # print(Pruned_T.vs[e.target]["name"])

                if v["penalties"] < e["cost"]:
                    # print("removo")
                    vertices_to_remove.append(v)

        if len(vertices_to_remove) == 0:
            break
        else:
            # removing vertices from vertices_to_remove and all its edges
            Pruned_T.delete_vertices(vertices_to_remove)

    return Pruned_T


def get_solution_list(E, T):
    print(E)
    S = [0] * len(E.keys())
    print(S)

    for e in T.es:
        link = (T.vs[e.source]["name"],T.vs[e.target]["name"])
        print("link")
        print(link)

        print("values of E")
        for (k, v) in E.items():
            print(v)
            if link in v:
                print(k)
                exit()

    return S


def get_graph_of_terminals():
    return


def create_xlsx(path, arg_number, instance, G, initial_S, star_S,
                initial_pont, star_pont, temp_ini, temp_f, f, SA_max, process_time, points):

    folder_name = 'Results_' + str(arg_number)
    if folder_name not in os.listdir(path):
        os.mkdir(os.path.join(path, folder_name))

    path_results = os.path.join(path, folder_name)

    file_name = instance.replace('.', '_').replace('-', '_') + '_result_{}.xlsx'.format(arg_number)

    file_path = os.path.join(path_results, file_name)

    writer = pd.ExcelWriter(file_path, engine='xlsxwriter')

    # ------------------------------------------------------------ X
    # For G

    G_vertex = G.get_vertex_dataframe()
    G_edges = G.get_edge_dataframe()

    V_guide = G_vertex['name'].to_dict()  # dicionário com {index do vertice em aux_G: Nome do vertice em aux_G}
    # dataframe com arestas de aux_G
    G_edges['source'].replace(V_guide, inplace=True)  # Substitui o valor de source (que e index) para NOME
    G_edges['target'].replace(V_guide, inplace=True)  # Substitui o valor de target (que e index) para NOME

    G_vertex.to_excel(writer, sheet_name='Instance', index=False, startrow=0, startcol=0)
    G_edges.to_excel(writer, sheet_name='Instance', index=False, startrow=len(G_vertex.index) + 2)

    # ------------------------------------------------------------ X

    # ------------------------------------------------------------ X
    # For initial solution


    initial_S_vertex = initial_S.get_vertex_dataframe()
    initial_S_edges = initial_S.get_edge_dataframe()

    V_guide = initial_S_vertex['name'].to_dict()  # dicionário com {index do vertice em aux_G: Nome do vertice em aux_G}
    # dataframe com arestas de aux_G
    initial_S_edges['source'].replace(V_guide, inplace=True)  # Substitui o valor de source (que e index) para NOME
    initial_S_edges['target'].replace(V_guide, inplace=True)  # Substitui o valor de target (que e index) para NOME

    initial_S_vertex.to_excel(writer, sheet_name='Initial_S', index=False, startrow=0, startcol=0)
    initial_S_edges.to_excel(writer, sheet_name='Initial_S', index=False, startrow=len(initial_S_vertex.index) + 2)

    # ------------------------------------------------------------ X
    # For star solution

    star_S_vertex = star_S.get_vertex_dataframe()
    star_S_edges = star_S.get_edge_dataframe()

    V_guide = star_S_vertex['name'].to_dict()  # dicionário com {index do vertice em aux_G: Nome do vertice em aux_G}
    # dataframe com arestas de aux_G
    star_S_edges['source'].replace(V_guide, inplace=True)  # Substitui o valor de source (que e index) para NOME
    star_S_edges['target'].replace(V_guide, inplace=True)  # Substitui o valor de target (que e index) para NOME

    star_S_vertex.to_excel(writer, sheet_name='Star_S', index=False, startrow=0, startcol=0)
    star_S_edges.to_excel(writer, sheet_name='Star_S', index=False, startrow=len(star_S_vertex.index) + 2)

    # ------------------------------------------------------------ X

    resume = pd.DataFrame({'Nome da Instancia': [instance], 'Pontuação Inicial': [initial_pont],
                           'Pontuação Final': [star_pont], 'T inicial': [temp_ini], 'T final': [temp_f],
                           'f decaimento': [f], 'iterações Metropolis': [SA_max],
                           'Tempo de Processamento': [process_time], 'Bilhetes': [points]})

    resume.to_excel(writer, sheet_name='Resume', index=False)

    writer.save()

    return path_results


def create_resume(results_path):

    files = os.listdir(results_path)

    df = pd.DataFrame()

    for f in files:

        f_path = os.path.join(results_path, f)

        df_aux = pd.read_excel(f_path, sheet_name='Resume')

        df = df.append(df_aux)

    out_dir = os.path.join(results_path, 'Resume.xlsx')
    df.to_excel(out_dir, index=False)

    return 'fim'


if __name__ == "__main__":
    dir = "instances/our-instances/"
    subdirs = os.listdir(dir)
    for subd in subdirs:

        path = dir + subd
        files = os.listdir(path)
        for f in files:
            filename = path + "/" + f
            print(filename)
            if subd == "pcstp":
                continue
            elif subd == "E" or subd == "SteinCD":
                G, V, vertex_penalties, E, edge_costs = create_graph_from(filename)
            else:
                G, V, vertex_penalties, E, edge_costs = create_graph_from_stp(filename)
            # print(V)
            # print(vertex_penalties)
            # print(E)
            # print(edge_costs)

            layout = G.layout("lgl")
            G.vs["label"] = G.vs["name"]
            plot(G,layout=layout)


    # arq = "instances/our-instances/PCSPG-ACTMODPC/drosophila005.stp"
    # G, V, vertex_penalties, E, edge_costs = create_graph_from_stp(arq)
    # print(V)
    # print(vertex_penalties)
    # print(E)
    # print(edge_costs)
    #
    # layout = G.layout("lgl")
    # G.vs["label"] = G.vs["name"]
    # plot(G,layout=layout)
