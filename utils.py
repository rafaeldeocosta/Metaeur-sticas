from igraph import Graph
from igraph import plot
import itertools as it


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
        E[int(e[0])] = ( int(e[1]), int(e[2]))
        edge_costs[int(e[0])] = int(e[3])


    G = Graph()
    G.add_vertices(len(V))  # adding the number of vertices of the graph
    G.vs["name"] = list(V) # setting labels to identify vertices

    # inserting edges of the graph
    for e in E.values():
        u = e[0]
        v = e[1]
        G.add_edge(G.vs.find(name=u), G.vs.find(name=v))

    G.es['cost'] = list(edge_costs.values())
    G.vs['penalties'] = list(vertex_penalties.values())

    return G, V, vertex_penalties, E, edge_costs


def select_sub_graph(G, vertices):

    """
    Função retorna um sub grafo com os mínimos caminhos entre os vértices da lista "vertices"
    :param G: igraph.Graph
        Grafo principal, onde o custo das arestas está no atributo nomeado de "cost" e as penalidades dos véritces em
        "penalties"
    :param vertices: list
        Lista de NOMES de vertices
    :return:
        G: igraph.Graph
            Subgrafo de G com vértices da lista "vertices"
        E_aux: dict
            Dicionario com {index da arestas: (u, v)} > u, v = nomes dos vértices source e target
        edges_cost_aux: dict
            Dicionario com {index arestas: cost aresta}
    """

    combi_V = list(it.combinations(vertices, 2))  # combinacao 2 a 2 dos itens na lista de vertices

    verts_G = G.get_vertex_dataframe().copy()

    aux_G = Graph()
    aux_G.vs['name'] = []

    for item in combi_V:

        source = item[0]
        destination = item[1]

        minimal_path = G.get_shortest_paths(G.vs.find(name=source), G.vs.find(name=destination), weights='cost')[0]  #retorna os indexes dos vértices do menor caminho

        for id, vert in enumerate(minimal_path[:-1]):

            G_u_idx = vert  # Index de u em G
            u = G.vs.find(vert)['name']  # Nome de u

            G_v_idx = minimal_path[id+1]  # Index de v em G
            v = G.vs.find(G_v_idx)['name']  # Nome de v

            if u not in aux_G.vs['name']:

                # aux_G.add_vertex(u)
                penalty_u = verts_G.loc[verts_G['name'] == u, 'penalties'].values[0]  # penalidade de u
                aux_G.add_vertices([u], attributes={'penalties': [penalty_u]})  # dessa maneira a penalidade já entra no G

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

if __name__ == "__main__":
    arq = "K100.1"
    G, V, vertex_penalties, E, edge_costs = create_graph_from(arq)
    print(V)
    print(vertex_penalties)
    print(E)
    print(edge_costs)

    layout = G.layout("lgl")
    G.vs["label"] = G.vs["name"]
    plot(G,layout=layout)
