#Dirichlet_CI
import math
from scipy.stats import multinomial
from decimal import Decimal
from PyQt6 import uic
from PyQt6.QtWidgets import QApplication, QWidget, QMainWindow
import plotly.express as px
import numpy as np
# from sympy.stats.joint_rv_types import MultivariateBeta
# from sympy.stats import density
# from sympy import Symbol, pprint





#param = [p_det, p_prom, p_pass]
#observation = [n_det, n_prom, n_pass]

def posterior_param_distribution(dirichlet_param, observation):
    # param is a list of the parameters of the Dirichlet distribution
    # observation is a list of the observations of the multinomial distribution
    # returns a list of the posterior parameters of the Dirichlet distribution
    posterior_param = []
    for i in range(0, len(dirichlet_param)):
        posterior_param.append(dirichlet_param[i] + observation[i])
    return posterior_param
    
def Beta(alpha_array):
    # alpha_array is a list of the parameters of the Dirichlet distribution
    # returns the Beta distribution corresponding to the Dirichlet distribution
    Beta = 1
    for i in range(0, len(alpha_array)):
        Beta = Beta * math.gamma(alpha_array[i])
    Beta = Beta / math.gamma(sum(alpha_array))
    return Beta

def Dirichlet_posterior_predictive(x, observation, dirichlet_param):
    # x is the value at which the posterior predictive distribution is evaluated
    # observation is a list of the observations of the multinomial distribution
    # dirichlet_param is a list of the parameters of the Dirichlet distribution
    # returns the value of the posterior predictive distribution at x
    #N = sum( x[i] for i in range(0, len(x)))
    N = 0
    for x_i in x:
        N = N + x_i
    #print(N)
    posterior_param = posterior_param_distribution(dirichlet_param, observation)
    posterior_param = posterior_param_distribution( posterior_param, x)
    posterior_predictive = math.factorial(N)*Beta(posterior_param)/Beta(dirichlet_param)
    for i in range(0, len(x)):
        posterior_predictive = posterior_predictive / math.factorial(x[i])
    
    return posterior_predictive

def proba_nps_knowing_observation(nps, N, observation, dirichlet_param):
    # probability of observing nps = y knowing the observation and given the fact that nps = x0 - x1 / (x0 + x1 + x2)
    # ie x0 = promoter
    # x1 = detractor
    # x2 = passive
    sum_proba = 0
    if nps >= 0:
        for i in range (0, int((N - nps)/2)+1):
            #print('Promoters {}, Detractors {}, Passive {}'.format(nps +i, i, N - nps - 2*i))
            x = [nps +i, i, N - nps - 2*i]
            sum_proba = sum_proba + Dirichlet_posterior_predictive(x, observation, dirichlet_param)
    
    #negative NPS case
    elif nps < 0:
        for i in range(0, int((N +nps)/2)+1):
            x = [i, -nps+i, N + nps - 2*i]  
            #print('Promoters {}, Detractors {}, Passive {}'.format(x[0], x[1], x[2]))
            sum_proba = sum_proba + Dirichlet_posterior_predictive(x, observation, dirichlet_param)
    

    else: return -1
    return sum_proba

def calculate_NPS(n_prom, n_det, N):
    nps = (n_prom - n_det)/N
    x = n_prom/N
    y = n_det/N
    return nps,x,y




def calculate_CI_augmented(n_prom, n_det, N, alpha, density_array):
    nps,x,y = calculate_NPS(n_prom, n_det, N)
    nps = round(nps, 2)
    #find which index of the NPS probabiity density records the observed nps
    index = - 100000000000000
    for i in density_array:
        if i[0]== nps:
            index = density_array.index(i)

    print ("index is {}".format(index))
    total_proba=density_array[index][1]
    for j in range (1, 2*N+1):
        if index+j < len(density_array):
            high_gift = density_array[index+j][1]
            high_bound_nps = density_array[index+j-1][0]
        else: high_gift=0

        if index - j >= 0:
            low_gift = density_array[index-j][1]
            low_bound_nps = density_array[index-j+1][0]
        else: low_gift=0

        if(total_proba + high_gift + low_gift >=1-alpha):

            
            return (high_bound_nps, low_bound_nps)
        total_proba = total_proba + high_gift + low_gift
        print('high {}, low {}, cumulated proba {}'.format(index+j, index-j, total_proba))
        print ("high bound {}, low bound {}".format(high_bound_nps,low_bound_nps))



if __name__ == "__main__":
    dirichlet_param = [1,1,1]
    observation = [10,2,3]
    x = [8,2,2]
    N = sum(observation[i] for i in range (0, len(observation)))
    print(N)
    print(Dirichlet_posterior_predictive(x, observation, dirichlet_param))
    

    total_sum = 0
    
    plot_y = np.array([])
    plot_x = np.array([])
    nps_distrib=[]
    for j in range (-12,13):
        result =proba_nps_knowing_observation(j, N, observation, dirichlet_param)
        total_sum = total_sum + result
        #print("NPS {} , cummulated probability {}, density {}".format(j, total_sum, result))
    norm_factor = total_sum
    #create an empty two dimensional np array
  

    for i in range (-12,13):
        result =proba_nps_knowing_observation(i, N, observation, dirichlet_param)
        total_sum = total_sum + result /norm_factor
        #plot_x = plot_x.append(round(i/12,2))
        plot_x = np.append(plot_x, round(i/N,2))
        plot_y = np.append(plot_y,round(result / norm_factor,5))
        nps_distrib.append((round(i/N,2),round(result / norm_factor,5)))
        print("NPS {} , cummulated probability {}, density {}".format(round(i/N,2), round(total_sum,3), round(result / norm_factor,5)))

    fig = px.bar(x=plot_x, y=plot_y)
    fig.show()
    max_likelihood_nps = np.argmax(plot_y)
    print ("Max likelihood NPS is {}".format(max_likelihood_nps))

    calculate_CI_augmented(10,2,15,0.10, nps_distrib)
