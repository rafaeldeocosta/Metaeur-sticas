from metropolis import metropolis
from utils import calc_pontuacao

def SA(G, Temp_ini, Temp_fin, ALPHA, S_ini, SA_MAX):
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

    # print(Temp_ini)
    # print(Temp_fin)
    # print(ALPHA)
    # print(S_ini)
    # print(SA_MAX)

    Temp_curr = Temp_ini
    Star_S = S_ini

    while(Temp_curr > Temp_fin):

        New_S = metropolis(G, S_ini, Temp_curr, SA_MAX)

        if calc_pontuacao(G, New_S) < calc_pontuacao(G, Star_S):
            Star_S = New_S


        Temp_curr = Temp_curr * ALPHA

    return Star_S
