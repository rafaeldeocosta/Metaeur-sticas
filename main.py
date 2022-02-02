from datetime import datetime
import os

from kruskal import kruskal
from SA import SA
from utils import calc_pontuacao
from utils import create_graph_from
from utils import remove_costly_leafs
from utils import select_sub_graph
from utils import vert_premios
from utils import create_xlsx
from utils import create_resume


if __name__ == "__main__":

    print('#-----------------------------------------------------------------#\n')
    print("Trabalho de Meta Heuristica: Simulated Annealing\n")
    print('#-----------------------------------------------------------------#\n')

    # Path of instances to get solution
    instances_path = 'intances_teste'

    # list of instances
    instances = os.listdir(instances_path)


    # TODO: get params from file called 'cooling_strategies'
    # [Temp_ini, Temp_fin, SA_max, f]
    args = [[1000, 1, 20, 'Temp_curr*0.5']]

    for arg in args:

        Temp_ini = arg[0]
        Temp_fin = arg[1]
        SA_MAX = arg[2]
        f = arg[3]  # Função decaimento da temperatura > IMPORTANTE: em função de "Temp_curr"

        arg_number = args.index(arg) + 1

        print('\n\n\n#-----------------------------------------------------------------#')
        print('#-----------------------------------------------------------------#')
        print('#-----------------------------------------------------------------#')
        print('#-----------------------------------------------------------------#\n')
        print("Set {} of arguments:\n Initial Temperature: {}\n"
              "Final Temperature: {}\n"
              "Number of iterations in Metropolis: {}\n"
              "Decay Function: {}".format(arg_number, Temp_ini, Temp_fin, SA_MAX, f))
        print('#-----------------------------------------------------------------#\n')

        for arq in instances:

            if 'Result' not in arq:

                print('#-----------------------------------------------------------------#\n')
                print("Instance being used: {}\n".format(arq))
                print('#-----------------------------------------------------------------#\n')


                # ###################################################################
                # Generate graph G#
                #####################################################################
                arq_path = os.path.join(instances_path, arq)

                G, V, vertex_penalties, E, edge_costs = create_graph_from(arq_path)

                # ###################################################################
                # Generate subgraph of G connecting terminals #
                #####################################################################

                terminais = vert_premios(G)

                G_terminais, E_terminais, edge_costs_terminais = \
                    select_sub_graph(G, terminais)

                # ##################################################
                # Generate partial initial solution with kruskal algorithm #
                ####################################################
                #

                T_terminais, T_terminais_edges_list = \
                    kruskal(G_terminais, E_terminais, edge_costs_terminais)

                sol_T_terminais = calc_pontuacao(G, T_terminais)

                print('#-----------------------------------------------------------------#\n')
                print("Score of Kruskal partial inicial solution: %s" % sol_T_terminais)

                # ##################################################
                # Generating initial solution removing costly leafs #
                ####################################################
                #

                Pruned_T_Terminais = remove_costly_leafs(T_terminais)

                sol_Pruned_T_terminais = calc_pontuacao(G, Pruned_T_Terminais)


                print("Score of initial solution: %s" % sol_Pruned_T_terminais)
                print('#-----------------------------------------------------------------#\n')


                # ##################################################
                # Running simulated anealing #
                ####################################################


                S = Pruned_T_Terminais.copy()

                start_time = datetime.now()

                Star_S = SA(G, Temp_ini, Temp_fin, S, SA_MAX, f)

                end_time = datetime.now()

                process_time = str(end_time - start_time)

                print('#-----------------------------------------------------------------#\n')
                print('Processing time: {}'.format(process_time))
                print('#-----------------------------------------------------------------#\n')

                sol_Star_S = calc_pontuacao(G, Star_S)

                print('#-----------------------------------------------------------------#\n')
                print("Score of Star_S %s" % sol_Star_S)
                print('#-----------------------------------------------------------------#\n')

                # ##################################################
                # Generating files with results #
                ####################################################

                path_results = create_xlsx(instances_path, arg_number, arq, G, S, Star_S, sol_Pruned_T_terminais, sol_Star_S,
                                        Temp_ini, Temp_fin, f, SA_MAX, process_time)


        create_resume(path_results)


    print('#-----------------------------------------------------------------#\n')
    print("FIM DO PROCESSAMENTO")
    print('#-----------------------------------------------------------------#\n')

