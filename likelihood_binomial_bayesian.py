#likelihood_binomial_bayesian.py
import numpy as np
import math 
import sys
from scipy.stats import multinomial
from decimal import Decimal
from PyQt6 import uic
from PyQt6.QtWidgets import QApplication, QWidget, QMainWindow


alpha_level = 0.15
resolution = 10


class UI(QMainWindow):

    def __init__(self):
        super().__init__()
        
        # loading the ui file with uic module
        uic.loadUi("C:\Repos Surface\Dev Perso\Stats\MainUI.ui", self)
        #self.label.setObjectName("Cuir")
        strVersion = 'Version 0'
        strComment = "The distribution of detractors and promoters is assumed to be multinomial. The probability for a detractor or promoters are calculated using the maximum likelihood method "
        self.label_version_comment.setText(strVersion + "-" + strComment)
        self.label_version.setText(strVersion)
        self.pushButton.clicked.connect(self.buttonclicked)
        self.progressBar.setHidden(True) 
        self.show()
    def buttonclicked(self):
        self.progressBar.setHidden(False)
        self.progressBar.setValue(0)

        arr1 = [0,0,0]
        arr2 = [0,0,0]
        arr3 = [0,0,0]
       
        
        try:
            arr1[0] = int(self.in_promoters_v1.toPlainText())
            arr1[1] = int(self.in_detractors_v1.toPlainText())
            arr1[2] = int(self.in_passives_v1.toPlainText())+int(self.in_detractors_v1.toPlainText())+int(self.in_promoters_v1.toPlainText())

            arr2[0] = int(self.in_promoters_v2.toPlainText())
            arr2[1] = int(self.in_detractors_v2.toPlainText())
            arr2[2] = int(self.in_passives_v2.toPlainText())+int(self.in_detractors_v2.toPlainText())+int(self.in_promoters_v2.toPlainText())

            print(self.in_promoters_v3.toPlainText())
            if self.in_promoters_v3.toPlainText() == "" and self.in_detractors_v3.toPlainText() == "" and self.in_passives_v3.toPlainText() == ""  :
                print("empty")
                arr3 = None

            else: 
                arr3[0] = int(self.in_promoters_v3.toPlainText())
                arr3[1] = int(self.in_detractors_v3.toPlainText())
                arr3[2] = int(self.in_passives_v3.toPlainText())+int(self.in_detractors_v3.toPlainText())+int(self.in_promoters_v3.toPlainText())
        
            testtry = arr1[0]*arr1[1]*arr1[2]
        except:
            print("Error, enter numbers please, not your twitter feed")
        print(arr1[0])
        print(arr1[1])
        print(arr1[2])

        if arr3 ==None:
            p1,p2,nps1, nps2 = do_the_test(arr1, arr2, arr3)
            self.nps_v1.setText(str(round(nps1,2)))
            self.nps_v2.setText(str(round(nps2,2)))

            self.pvalue_v1.setText(str(round(p1,2)))
            self.pvalue_v2.setText(str(round(p2,2)))

            if p1 < alpha_level or p2 < alpha_level:
                
                self.label_result.setStyleSheet('QWidget{\ncolor: rgb(10, 255, 10);}')
                self.label_result.setText("There is a significant difference between the populations at 0.85 confidence level")
            else: 
                self.label_result.setText("There is no significant difference between populations at 0.85 confidence level ")
        else:
            p1,p2,p3,nps1, nps2, nps3 = do_the_test(arr1, arr2, arr3)
            self.nps_v1.setText(str(round(nps1,2)))
            self.nps_v2.setText(str(round(nps2,2)))
            self.nps_v3.setText(str(round(nps3,2)))

            self.pvalue_v1.setText(str(round(p1,2)))
            self.pvalue_v2.setText(str(round(p2,2)))
            self.pvalue_v3.setText(str(round(p3,2)))

            if p1 < alpha_level or p2 < alpha_level or p2 < alpha_level:
                    
                self.label_result.setStyleSheet('QWidget{\ncolor: rgb(10, 255, 10);}')
                self.label_result.setText("There is a significant difference between the populations at 0.85 confidence level")
            else: 
                self.label_result.setText("There is no significant difference between populations at 0.85 confidence level ")
        


def NPS_probability(n_prom, n_det, N, x,y):
    sum = 0
    nps = n_prom - n_det
    for i in range (0,N-n_prom+n_det+1):
        nps = n_prom-(n_det)+i
        for j in range (0,int((N-nps)/2)+1):
            sum = sum + mymultinomial(nps+j,j,N,x,y)
    return(sum)

def mymultinomial(n_prom,n_det,N,x,y):
    # calculates the probability of having a given numbers of promoters, detractors and passives in an NPS
    # survey with N responders, given the probability for a responder of being a promoter, y of being a detractor
    # and 1-x-y of being passive. 
    #print("Compute probability of multibinomial distribution of the following params: ")
    if (N-n_prom-n_det) <0 or 1-x-y<0:
        return 0
    else:

        result = multinomial(N,[x,y,1-x-y]).pmf([n_prom, n_det, N-n_prom-n_det])
        return result

def cumulated_proba(n_prom, n_det, N, x,y):
    sum = 0
    for i in range(0, N):
        for j in range(0,N):
            sum = sum+mymultinomial(i,j,N,x,y)
    print("Cumulated proba computed: ")
    return(sum)

def Likelihood(x,y, arr1, arr2, arr3=None):
    one = mymultinomial(arr1[0],arr1[1],arr1[2],x,y)
    two = mymultinomial(arr2[0],arr2[1],arr2[2],x,y)
    if arr3 == None: 
        three = 1
    else: 
        three = mymultinomial(arr3[0],arr3[1],arr3[2],x,y)

    
    return one*two*three


#def MaxLikelihood(a1,b1,N1,a2,b2,N2,a3,b3,N3):
def MaxLikelihood(arr1, arr2, arr3=None):



    maxL = 0
    best_x =0
    best_y = 0
    ultimate_x=0
    ultimate_y=0
    for i in range (0,100):
        window.progressBar.setValue(i/2) 
        for j in range (0,100):
            result = Likelihood(i/100,j/100,arr1,arr2,arr3)
            if result > maxL:
                maxL = result
                best_x = i/100 
                best_y = j/100
    ultimate_x=best_x
    ultimate_y=best_y
    for i in range(-100,100):
        window.progressBar.setValue(50+(100+i)/4) 
        for j in range (-100,100):
            result = Likelihood(best_x+i/1000,best_y+j/1000,arr1,arr2,arr3)
            if result > maxL:
                maxL = result
                ultimate_x = best_x+i/1000
                ultimate_y = best_y+i/1000

    return (maxL, ultimate_x, ultimate_y)

def probability_of_observation_given_xy(x,y,arr1, arr2, arr3=None, arr4=None):
    # computes the probability of having a given observation, given x and y, probabilities of being a promoter and detractor
    # arr1, arr2, arr3, arr4 are arrays of the form [n_promoters, n_detractors, n_passives]
    # arr4 is optional, if not provided, the function will compute the probability of arr1 and arr2 only
    # the function returns the probability of having the given observation, given the probabilities of being a promoter or detractor
    # the function will also return the probability of having the given observation, given the probabilities of being a promoter or detractor, for all possible combinations of x,y, probabilities of being a 
    # promoter or detractor

    probability =0 
    factor1 = mymultinomial(arr1[0],arr1[1],arr1[2],x,y)
    factor2 = mymultinomial(arr2[0],arr2[1],arr2[2],x,y)
    if arr3 == None:
        factor3 = 1
    else:
        factor3 = mymultinomial(arr3[0],arr3[1],arr3[2],x,y)
    if arr4 == None:
        factor4 = 1 
    else:
        factor4 = mymultinomial(arr4[0],arr4[1],arr4[2],x,y)
    probability = factor1*factor2*factor3*factor4
    
    return probability

def probabilty_of_observation_given_binomial_distribution(arr1, arr2, arr3=None, arr4=None):
    #it is the sum of probabilities of having the given observation, given the probabilities of being a promoter or detractor, for all possible combinations of x,y
    sum = 0
    for k, h in np.ndindex((resolution,resolution)):
        pk=k/resolution
        ph=h/resolution
        factor1 = mymultinomial(arr1[0],arr1[1],arr1[2],pk,ph)
        factor2 = mymultinomial(arr2[0],arr2[1],arr2[2],pk,ph)
        if arr3 == None:
            factor3 = 1
        else:
            factor3 = mymultinomial(arr3[0],arr3[1],arr3[2],pk,ph)
        if arr4 == None:
            factor4 = 1
        else:
            factor4 = mymultinomial(arr4[0],arr4[1],arr4[2],pk,ph)
        probability = factor1*factor2*factor3*factor4
        sum = sum + probability
    return sum




def probability_of_p_d_p(n_prom, n_det, N,arr1,arr2,arr3=None,arr4=None):
    # computes the probability of having a given number of promoters and detractors, given a total number N of responders, over all possible combinations of x,y, probabilities of being a 
    # promoter or detractor, and given observations arr1, arr2, arr3 and arr4

    sum = 0
    for k, h in np.ndindex((resolution,resolution)):
        pk=k/resolution
        ph=h/resolution
        
        sum = sum + mymultinomial(n_prom,n_det,N,pk,ph) * probability_of_observation_given_xy(pk,ph,arr1,arr2,arr3,arr4)
    #window.progressBar.setValue(k)
    
    result = sum / probabilty_of_observation_given_binomial_distribution(arr1,arr2,arr3,arr4)
    #print ("probability of having {} promoters and {} detractors, given {} total responders, is {}".format(n_prom,n_det,N,result))
    return result


def proba_nps(NPS,N,arr1,arr2,arr3 = None, arr4=None):
    sum = 0
    for i in range(0,int((N - NPS)/2)+1):
        sum = sum + probability_of_p_d_p(NPS+i, i,N,arr1,arr2,arr3, arr4)
    print ("Probability of nps: {} is {}".format(NPS,sum))
    return sum       

                
def do_the_test(arr1,arr2,arr3 = None):
        (p, x, y) = MaxLikelihood(arr1,arr2,arr3)       
        print("Probability: {}, p_promoter: {}, p_detractor: {}".format(p,x,y))
    
        p_value_v1 = NPS_probability(arr1[0],arr1[1],arr1[2],x,y)
        if p_value_v1 >0.5:  p_value_v1 = 1 - p_value_v1
        print ("p-value for first variant: {}, NPS first variant: {} ".format(p_value_v1,(arr1[0] - arr1[1])/arr1[2]))
        

        p_value_v2 = NPS_probability(arr2[0],arr2[1],arr2[2],x,y)
        window.progressBar.setHidden(True)
        if p_value_v2 >0.5: p_value_v2 = 1 - p_value_v2
        print ("p-value for second variant: {}, NPS second variant: {}".format(p_value_v2,(arr2[0] - arr2[1])/arr2[2]))
        
    

        if arr3 == None:
            print("Only two variants")
            return p_value_v1, p_value_v2, (arr1[0] - arr1[1])/arr1[2], (arr2[0] - arr2[1])/arr2[2]
        else:
            p_value_v3 = NPS_probability(arr3[0],arr3[1],arr3[2],x,y)
            if p_value_v3 >0.5: p_value_v3 = 1 - p_value_v3
            print ("p-value for third variant: {}, NPS third variant: {}".format(p_value_v3,(arr3[0] - arr3[1])/arr3[2]))
            return p_value_v1, p_value_v2, p_value_v3, (arr1[0] - arr1[1])/arr1[2], (arr2[0] - arr2[1])/arr2[2], (arr3[0] - arr3[1])/arr3[2]
         




            
                
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = UI()
    # app.exec()
    arr1= [17,4,37]
    arr2= [20,6,35]
    N=6
    NPS = 4
    sum = 0
    for j in range(-N,N+1):
        sum = sum + proba_nps(j,N,arr1,arr2)
    print("Total sum of probabilities is {}, should be 1".format(sum))


    #do_the_test([17,4,37],[20,6,35])

   
