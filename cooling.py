from math import log

def get_cooling_strategy(cooling_str):
    """
        function that returns f equation to

        https://www.cos.ufrj.br/~daniel/mcmc/slides/aula_15.pdf (pag. 14)
    """

    if cooling_str == "exponential":
        f = 'Temp_ini*pow(0.5,t)'
    elif cooling_str == "linear":
        f = 'Temp_ini - (10 * t)'
    elif cooling_str == "logarithm":
        f = '1 / log(t + 0.5)'
    elif cooling_str == 'geometric':
        f = 'Temp_curr*0.9'

    return f

if __name__ == "__main__":
    f = get_cooling_strategy("logarithm")
    Temp_ini = 100
    t = 1
    ans = eval(f)
    print(ans)
