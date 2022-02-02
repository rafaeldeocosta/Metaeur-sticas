from cooling import get_cooling_strategy
from metropolis import metropolis
from utils import calc_pontuacao


def SA(G, Temp_ini, Temp_fin, S_ini, SA_MAX, cooling_str, type='SA'):
    """
        função que implementa o simulated annealing

        Args
            G - igraph.Graph - original Graph
            S_ini - igraph.Graph - Solution of PCSTP ("Steiner" Tree)
            Temp_ini - float - Initial Temperature
            Temp_fin - float - Final Temperature
            ALPHA - float - Cooling Rate
            SA_MAX - int - number of  iteration in the metropolis function
        Returns

    """
    print("Starting SA ...")

    # Number of steps
    t = 0

    # print(Temp_ini)
    # print(Temp_fin)
    # print(ALPHA)
    # print(S_ini)
    # print(SA_MAX)

    Temp_curr = Temp_ini
    Star_S = S_ini.copy()
    New_S = S_ini.copy()

    points = [1,1,1]

    if type == 'SA':
        incr = 0
    if type == 'SA-LNS':
        incr = 10

    while(Temp_curr > Temp_fin):

        print(points)

        Best_S, points, Best_process = metropolis(G, New_S,
                                                    Temp_curr, Temp_ini,
                                                    SA_MAX, points, type)


        if calc_pontuacao(G, Best_S) < calc_pontuacao(G, Star_S):
            Star_S = Best_S.copy()
            if Best_process == 'remove_vertices':
                points[0]+=incr
            elif Best_process == 'remove_edges':
                points[1]+=incr
            elif Best_process == 'both':
                points[2]+=incr

        New_S = Best_S.copy()

        # Cooling strategy
        f = get_cooling_strategy(cooling_str)
        Temp_curr = eval(f)

        t+=1

    return Star_S, points
