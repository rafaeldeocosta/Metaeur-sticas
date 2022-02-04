from utils import join_resumes
from utils import graph_results


if __name__ == "__main__":

    join_resumos = True
    graficos = False

    # Path é a pasta onde estão os grupos rodados, vale ressaltar que dentro dela devem ter PASTAS de grupos, os nomes
    # dessas pastas vão ser usadas para compor o resumo final
    if join_resumos:
        path = 'instances/rodadas_wagner'
        join_resumes(path)

    # Pasta com conjunto de instâncias que se deseja gerar gráficos, não coloque
    if graficos:
        path_instances = 'instances/rodadas_wagner/escolhidos'
        graph_results(path_instances)