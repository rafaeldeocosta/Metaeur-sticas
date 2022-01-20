from metropolis import metropolis

def SA(Temp_ini, Temp_fin, ALPHA, S_ini, SA_MAX):
    """
        função que implementa o simulated annealing

        Args:

            S_ini - igraph.Graph - Solution of PCSTP ("Steiner" Tree)
            Temp_ini - float - Initial Temperature
            Temp_fin - float - Final Temperature
            ALPHA - float - Cooling Rate
            SA_MAX - int - number of  iteration in the metropolis function
        Returns:

    """
    Temp_curr = Temp_ini
    Best_S = S_ini

    while(Temp_curr > Temp_fin):

        metropolis(S_ini, Temperature, SA_MAX)

        Temp_curr = Temp_curr * ALPHA




    return Best_S
