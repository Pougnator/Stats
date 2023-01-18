import math
from scipy.stats import multinomial
from decimal import Decimal
from PyQt6 import uic
from PyQt6.QtWidgets import QApplication, QWidget, QMainWindow
#Form,Window = uic.loadUiType("MainUI.ui")
# app = QApplication([])
# window = Window()
# form = Form()
#form.setupUi(window)



class UI(QMainWindow):

    def __init__(self):
        super().__init__()
        
        # loading the ui file with uic module
        uic.loadUi("MainUI.ui", self)
        #self.label.setObjectName("Cuir")
        self.pushButton.clicked.connect(self.buttonclicked)
    def buttonclicked(self):
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
        else:
            p1,p2,p3,nps1, nps2, nps3 = do_the_test(arr1, arr2, arr3)
            self.nps_v1.setText(str(round(nps1,2)))
            self.nps_v2.setText(str(round(nps2,2)))
            self.nps_v3.setText(str(round(nps3,2)))

            self.pvalue_v1.setText(str(round(p1,2)))
            self.pvalue_v2.setText(str(round(p2,2)))
            self.pvalue_v3.setText(str(round(p3,2)))


        


def NPS_probability(n_prom, n_det, N, x,y):
    sum = 0
    nps = n_prom - n_det
    for i in range (0,N-n_prom+n_det+1):
        nps = n_prom-(n_det)+i
        #print(nps)
        #print(nps)
        for j in range (0,int((N-nps)/2)+1):
            sum = sum + mymultinomial(nps+j,j,N,x,y)
        #print(sum)
    #print(sum)
    return(sum)

def mymultinomial(n_prom,n_det,N,x,y):
    # calculates the probability of having a given numbers of promoters, detractors and passives in an NPS
    # survey with N responders, given the probability for a responder of being a promoter, y of being a detractor
    # and 1-x-y of being passive. 
    #print("Compute probability of multibinomial distribution of the following params: ")
    if (N-n_prom-n_det) <0 or 1-x-y<0:
        return 0
    else:
        #result2 = stats.multinomial(N,(x,y,1-x-y))
        result = multinomial(N,[x,y,1-x-y]).pmf([n_prom, n_det, N-n_prom-n_det])
        #result =pow(x,n_prom)*pow(y,n_det)*pow((1-x-y),N-n_prom-n_det)*math.factorial(N)/math.factorial(n_prom)/math.factorial(n_det)/math.factorial(N-n_prom-n_det)
        #print (result)
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
   # print(one)
   # print(two)
   # print(three)
    
    return one*two*three

#def MaxLikelihood(a1,b1,N1,a2,b2,N2,a3,b3,N3):
def MaxLikelihood(arr1, arr2, arr3=None):

    # a1 = arr1[0]
    # b1 = arr1[1]
    # N1 = arr1[2]

    # a2 = arr2[0]
    # b2 = arr2[1]
    # N2 = arr2[2]
    # if arr3 == None:

    # else: 

 

    # a3 = arr3[0]
    # b3 = arr3[1]
    # N3 = arr3[2]

    maxL = 0
    best_x =0
    best_y = 0
    ultimate_x=0
    ultimate_y=0
    for i in range (0,100):
        for j in range (0,100):
            result = Likelihood(i/100,j/100,arr1,arr2,arr3)
            if result > maxL:
                maxL = result
                best_x = i/100 
                best_y = j/100
    ultimate_x=best_x
    ultimate_y=best_y
    for i in range(-100,100):
        for j in range (-100,100):
            result = Likelihood(best_x+i/1000,best_y+j/1000,arr1,arr2,arr3)
            if result > maxL:
                maxL = result
                ultimate_x = best_x+i/1000
                ultimate_y = best_y+i/1000

    return (maxL, ultimate_x, ultimate_y)
                
def do_the_test(arr1,arr2,arr3 = None):
        (p, x, y) = MaxLikelihood(arr1,arr2,arr3)        
        print("Probability: {}, p_promoter: {}, p_detractor: {}".format(p,x,y))
    
        p_value_v1 = NPS_probability(arr1[0],arr1[1],arr1[2],x,y)
        if p_value_v1 >0.5:  p_value_v1 = 1 - p_value_v1
        print ("p-value for first variant: {}, NPS first variant: {} ".format(p_value_v1,(arr1[0] - arr1[1])/arr1[2]))
        

        p_value_v2 = NPS_probability(arr2[0],arr2[1],arr2[2],x,y)
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


app = QApplication([])
window = UI()
window.show()
app.exec()



            
                
        
if __name__ == "__main__":
    # print("Multibinomal result:")
    # print(multinomial (8,1,15,8/15,1/15))
    # print(cumulated_proba(8,1,15,8/15,1/15))
    #print("likelihood: ")
    #print(Likelihood(8,1,15,11,5,18,11,5,20,0.562,0.21))
    #print(MaxLikelihood(8,1,15,11,5,18,11,5,20))
    # (p, x, y) = MaxLikelihood(13,7,22,11,6,20,13,12,30)
    # print(p,x,y)
    # NPS_probability(13,7,22,x,y)
    # NPS_probability(11,6,20,x,y)
    # NPS_probability(13,12,30,x,y)
    # print((13-2)/22)
    # print((11-3)/20)
    # print((13-5)/30)
    
    #do_the_test([13,7,22],[11,6,20],[13,12,30])
    do_the_test([17,4,37],[20,6,35])

   
