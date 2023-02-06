#Calculate NPS confidence interval

import math
from scipy.stats import multinomial
from decimal import Decimal
from PyQt6 import uic
from PyQt6.QtWidgets import QApplication, QWidget, QMainWindow
import plotly.express as px
import numpy as np




def NPS_density(N, x,y):
    # returns the probability density for NPS scores, given N respondents and probability x and y
    # for a respondent of being a promoter or a detractor respectively
    # returns an array of NPS,probability couples
    
    npsxN = []
    plot_y = []
    plot_x = []
    for i in range (-N,N+1):
        #i is the current NPSxN
        sum_proba =0
        n_prom = i
        for j in range(0, int((N - i)/2)+1):
            #The highest number of detractors possible for a given NPSxN is 
            #max_prom - max_det = i
            # and max_prom + max_det = N (because no more passives left)
            # so i + jmax + jmax = N
            # and jmax = (N-i)/2 of N-i even and (N-i-1)/2 if N-i is event (one passive left in this case)
            # both of this cases are summarised by the formula int((N-i)/2)
            n_det = j
            n_prom = i+j
            proba = mymultinomial(n_prom,n_det,N,x,y)
            sum_proba=sum_proba+proba
        npsxN.append((i/N,sum_proba))
        plot_y.append(sum_proba)
        plot_x.append(i/N)
    fig = px.bar(x=plot_x, y=plot_y)
    fig.show()
    return npsxN


    # sum = 0
    # nps = n_prom - n_det
    # for i in range (0,N-n_prom+n_det+1):
    #     nps = n_prom-(n_det)+i
    #     #print(nps)
    #     #print(nps)
    #     for j in range (0,int((N-nps)/2)+1):
    #         sum = sum + mymultinomial(nps+j,j,N,x,y)
    #     #print(sum)
    # #print(sum)
    # return(sum)


def mymultinomial(n_prom,n_det,N,x,y):
    # returns the probability of having a given numbers of promoters, detractors and passives in an NPS
    # survey with N responders, given the probability for a responder of being a promoter, y of being a detractor
    # and 1-x-y of being passive. 
    
    if (N-n_prom-n_det) <0 or 1-x-y<0:
        return 0
    else:
        
        result = multinomial(N,[x,y,1-x-y]).pmf([n_prom, n_det, N-n_prom-n_det])
        return result


def calculate_NPS(n_prom, n_det, N):
    nps = (n_prom - n_det)/N
    x = n_prom/N
    y = n_det/N
    return nps,x,y

def calculate_CI(n_prom, n_det, N, alpha):
    nps,x,y = calculate_NPS(n_prom, n_det, N)
    #find which idex of the NPS probabiity density records the observed nps
    density_array = NPS_density(N, x, y)
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

def calculate_CI_augmented(n_prom, n_det, N, alpha, density_array):
    nps,x,y = calculate_NPS(n_prom, n_det, N)
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

def bayesian_caluc(n_prom, n_det, N, alpha):
    nps,xobs,yobs = calculate_NPS(n_prom, n_det, N)
    #Calculating the probability P[(x,y) = j,k | xobs, yobs] 
    #P[(x,y) = j,k | xobs, yobs] = P(xobs, yobs | (x,y)=j,k) *P(x,y)=j,k / P(xobs,yobs)
    #[P(x,y)=j,k] = 1 because we consider these follow a uniform distrubution a priori
    #P(xobs,yobs) = sum P(xobs,yobs | x,y = j,k) over all possible j,k
    #And P(xobs, yobs) | (x,y)=j,k = multinomial(prom, det, N, j,k)
    # 
    #P(nps |xobs,yobs) =sum(P(nps|x=k,y=j,xobs,yobs)*P[(x,y) = j,k | xobs, yobs]
    #So finally: 
    #P(nps |xobs,yobs) =sum(P(nps|x=k,y=j,xobs,yobs)* P(xobs,yobs|(x,y)=j,k)/sum P(xobs,yobs|x,y = l,m) over all (l,m)
def NPS_density_augmented(N, oprom,odet):
    # returns the probability density for NPS scores, given N respondents and probability x and y
    # for a respondent of being a promoter or a detractor respectively
    # returns an array of NPS,probability couples
    
    npsxN = []
    plot_y = []
    plot_x = []
    for i in range (-2,N+1):
        #i is the current NPSxN
        sum_proba =0
        n_prom = i
        bfactor=0
        for k, h in np.ndindex((100,100)):
            pk=k/100
            ph=h/100
            sum_proba_nps =0
            for j in range(0, int((N - i)/2)+1):
                #The highest number of detractors possible for a given NPSxN is 
                #max_prom - max_det = i
                # and max_prom + max_det = N (because no more passives left)
                # so i + jmax + jmax = N
                # and jmax = (N-i)/2 of N-i even and (N-i-1)/2 if N-i is event (one passive left in this case)
                # both of this cases are summarised by the formula int((N-i)/2)
                n_det = j
                n_prom = i+j
                proba = mymultinomial(n_prom,n_det,N,pk,ph)
                sum_proba_nps=sum_proba_nps+proba
            #multiply by bayesian factor 
            bfactor = bfactor+mymultinomial(oprom,odet,N,pk,ph)
            sum_proba=sum_proba_nps*mymultinomial(oprom,odet,N,pk,ph)+sum_proba
        print("NPS {}, B factor {}".format(i,bfactor))
        npsxN.append((i/N,sum_proba/bfactor))
        plot_y.append(sum_proba/bfactor)
        plot_x.append(i/N)
    fig = px.bar(x=plot_x, y=plot_y)
    fig.show()
    return npsxN


if __name__ == "__main__":
    test_array = NPS_density_augmented(12, 9, 1)
    print(test_array)
    msum = 0
    for x in test_array:
        msum = msum + x[1]
        print (msum)
    #calculate_CI_augmented(9,1,12,0.05, test_array)