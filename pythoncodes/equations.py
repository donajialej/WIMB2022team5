import numpy as np
#import scipy.linalg as la
from scipy.optimize import fsolve



def solve_difference_equations(p, betat):

    '''Solves difference equations for our model'''
    

    time_vector = range(p.maxtime)

    old_max_infected_age = p.max_infected_age
    if p.simple:
        p.max_infected_age = p.c

    # Create arrays for susceptible (S), recovered (R), and deceased (D) individuals
    S, R, D = np.zeros(p.maxtime + 1), np.zeros(p.maxtime + 1), np.zeros(p.maxtime + 1)

    # Create arrays for daily (Iday) and total (Itot) infecteds
    Iday = np.zeros((p.max_infected_age, p.maxtime + 1))  # Get shape right
    Itot = np.zeros(p.maxtime + 1)

    # Initial number of recovereds and deceased, Rinit and Dinit
    Rinit, Dinit = 0, 0
    
    # Initialise Iday, Itot
    if p.simple:
        for j in range(p.max_infected_age):
            # Approximate distribution as linear over 21 days
            Iday[j, 0] = Iinit1 - j * (Iinit1 / (p.max_infected_age - 1))
    else:
        pseudo_dist = pseudo_IC_dist(p, rho)
        Iday[:, 0] = Iinit1 * pseudo_dist
        #Iday[0, 0] = Iinit1
    Itot[:] = 0
    Iinit = int(np.sum(Iday[:, 0]))

    # Initial number of susceptible (everyone that isn't already infected/recovered/dead)
    Sinit = p.N - (Iinit + Rinit + Dinit)

    # Initialise S, R, D
    S[0], R[0], D[0] = Sinit, Rinit, Dinit

    if travel_data:
        if p.square_lockdown:
            alpha = step(p, lgoog_data=len(p.alpha), parameters_dictionary=parameters_dictionary)
        else:
            alpha = tanh_spline(p, lgoog_data=len(p.alpha), parameters_dictionary=parameters_dictionary)

    if p.calculate_rate_vectors:
        p.beta, p.gamma, p.zeta = make_rate_vectors(parameters_dictionary, p)
        # print('making rate vectors...')

    A = np.zeros((p.max_infected_age, p.max_infected_age))
    A[:, 0] = rho * p.beta[0,:]
    for i in range(p.max_infected_age - 1):
        A[ i, i + 1] = 1 - p.gamma[0, i] - p.zeta[0,i]


    for i in time_vector:
        if travel_data:
            A[:, 0] = alpha[i] * rho * p.beta[0,:]
        # Compute total number of infecteds
        Infecteds = Iday[:, i]
        Itot[i] = np.sum(Infecteds)
        Iday[:, i + 1] = Infecteds @ A
        Iday[0, i + 1] = (S[i] / p.N) * Iday[0, i + 1]
        S[i + 1] = S[i] -  Iday[0, i+1] #(S[i] / p.N) * pressure[i]
        R[i + 1] = R[i] + np.dot(p.gamma, Infecteds)
        # Compute death rate per day (number of deceased would be D[i+1] = D[i] + np.dot(zeta, Infecteds))
        D[i + 1] = np.dot(p.zeta, Infecteds)
    Itot[p.maxtime] = np.sum(Iday[:, p.maxtime])

    p.max_infected_age = old_max_infected_age
    
    return S, Iday, R, D, Itot


