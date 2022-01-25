from igraph import Graph
from math import inf
from random import randint

def dijkstra(G, source):
    """
        função que implementa o algoritmo de Dijkstra para calcular o caminho
            mais curto entre dois vértices

            https://www.cos.ufrj.br/~daniel/grafos/slides/aula_8.pdf

        Args:
            G - igraph.Graph - grafo
            source - int - nome do vertice

        Returns:
            dist - dict - tamanho dos caminhos minimos de source até os demais
                            vertices de G; a key é o vertice e value é o tamanho

    """
    V = []                    # list of vertices of G
    for v in G.vs:
        V.append(v["name"])


    dist = {}                   # distance of source to each key of dist dict
    for v in V:
        dist[v] = inf

    dist[source] = 0

    explored_vertices = []      # list of explored vertices

    while sorted(explored_vertices) != sorted(V):

        # Selecione  u  em  D = V - S,  tal  que  dist[u]  é  mínima
        D = list(set(V) - set(explored_vertices))

        if len(D) == 0:
            continue

        u = None
        for u_ in D:
            if u is None:
                u = u_
            elif dist[u_] < dist[u]:
                u = u_
            else:
                continue

        # Adicione  u  em  S
        explored_vertices.append(u)

        # Para  cada  vizinho  v  de  u  faça
        adj_vertices = G.neighbors(u)
        viz = G.vs[adj_vertices]["name"]
        for v in viz:
            w_u_v = G.es.select(_source=u, _target=v)["cost"][0]
            if dist[v] is inf or dist[v] > dist[u] + w_u_v:
                dist[v] = dist[u] + w_u_v

    print(explored_vertices)
    return dist

if __name__ == "__main__":

    """
        Exemplo retirado da aula 8 - Grafos com pesos, comprimento, distâncias,
            ideia e algoritmo de Dijkstra, Dijkstra - o próprio.
            https://www.cos.ufrj.br/~daniel/grafos/slides/aula_8.pdf
    """

    G = Graph()
    G.add_vertices(7)
    G.vs["name"] = ["a","b","c","d","e","f","g"]

    # A - B
    G.add_edge(G.vs.find(name="a"),G.vs.find(name="b"))
    e = G.es.select(_source="a", _target="b")
    e["cost"] = [1]

    # A - D
    G.add_edge(G.vs.find(name="a"),G.vs.find(name="d"))
    e = G.es.select(_source="a", _target="d")
    e["cost"] = [5]

    # A - E
    G.add_edge(G.vs.find(name="a"),G.vs.find(name="e"))
    e = G.es.select(_source="a", _target="e")
    e["cost"] = [2]

    # B - C
    # G.add_edge(G.vs.find(name="b"),G.vs.find(name="c"))
    # e = G.es.select(_source="b", _target="c")
    # e["cost"] = [4]

    # B - D
    G.add_edge(G.vs.find(name="b"),G.vs.find(name="d"))
    e = G.es.select(_source="b", _target="d")
    e["cost"] = [2]

    # C - D
    # G.add_edge(G.vs.find(name="c"),G.vs.find(name="d"))
    # e = G.es.select(_source="c", _target="d")
    # e["cost"] = [1]

    # C - G
    G.add_edge(G.vs.find(name="c"),G.vs.find(name="g"))
    e = G.es.select(_source="c", _target="g")
    e["cost"] = [2]

    # D - E
    G.add_edge(G.vs.find(name="d"),G.vs.find(name="e"))
    e = G.es.select(_source="d", _target="e")
    e["cost"] = [2]

    # D - F
    # G.add_edge(G.vs.find(name="d"),G.vs.find(name="f"))
    # e = G.es.select(_source="d", _target="f")
    # e["cost"] = [2]

    # E - F
    # G.add_edge(G.vs.find(name="e"),G.vs.find(name="f"))
    # e = G.es.select(_source="e", _target="f")
    # e["cost"] = [3]

    # F - G
    G.add_edge(G.vs.find(name="f"),G.vs.find(name="g"))
    e = G.es.select(_source="f", _target="g")
    e["cost"] = [3]

    # dist = dijkstra(G, "a")
    # print(dist)
    # c = G.clusters()
    # print(len(c))
    # s1 = c.subgraph(0)
    # print(type(s1))
    # print(c[0].subgraph())
    # print(c[1].subgraph())

    print(G.get_vertex_dataframe())
    u = G.vs[randint(0, len(G.vs) - 1)]
    v = G.vs[randint(0, len(G.vs) - 1)]
    print(u["name"])
    print(v["name"])

    P = G.get_all_shortest_paths(G.vs.find(name=u["name"]),
                            to=G.vs.find(name=v["name"]))

    print(P)
    print(len(G.es.select(_source=u["name"], _target=v["name"])))
    # print(e[0])
    # print(e[0].source)
    # print(e[0].target)
