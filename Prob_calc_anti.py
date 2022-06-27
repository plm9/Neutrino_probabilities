import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from tqdm import tqdm
import seaborn as sns

from Prob_calc import PMNS_param_matrix,flavor_to_index,c_ij,s_ij

#Defining symbols
s_12=Symbol("s_12",real=True)
s_23=Symbol("s_23",real=True)
s_13=Symbol("s_13",real=True)

c_12=Symbol("c_12",real=True)
c_23=Symbol("c_23",real=True)
c_13=Symbol("c_13",real=True)

d_cp=Symbol("d_cp",real=True)

D_m_12=Symbol("D_m_12",real=True)
D_m_13=Symbol("D_m_13",real=True)
D_m_23=Symbol("D_m_23",real=True)

L=Symbol("L",real=True,positive=True)
E_nu=Symbol("E_nu",real=True)

PMNS=PMNS_param_matrix()

def Prob_anti_a_to_anti_b(a,b):
    index_a,index_b=flavor_to_index(a,b)

    f=conjugate(PMNS[index_b,0])*PMNS[index_a,0]*conjugate(PMNS[index_a,1])*PMNS[index_b,1]
    s=conjugate(PMNS[index_b,0])*PMNS[index_a,0]*conjugate(PMNS[index_a,2])*PMNS[index_b,2]
    t=conjugate(PMNS[index_b,1])*PMNS[index_a,1]*conjugate(PMNS[index_a,2])*PMNS[index_b,2]

    first=re(f)*sin((D_m_12*L)/(4*E_nu))**2
    second=re(s)*sin((D_m_13*L)/(4*E_nu))**2
    third=re(t)*sin((D_m_23*L)/(4*E_nu))**2

    if a==b:
        return 1-4*(first+second+third)
    else:
        return -4*(first+second+third)

delta_cp=300

e="e"
mu="mu"
tau="tau"

flavors=[e,mu]

#These variables are used fot the plots
n=200
dist=np.linspace(0,13000,n)
eneg_range=np.logspace(-1,2,n)

list_of_all=[]

#Calculating the probabilities for all flavors
bool=True
for a in tqdm(flavors):
    for b in flavors:       
        if bool or a==b:
            #Replacing sin and cos by their values known from experiments
            Total_prob=np.zeros([n,n])
            Equation=Prob_anti_a_to_anti_b(a,b).subs([(s_12,s_ij(1,2)),(c_12,c_ij(1,2)),(c_13,c_ij(1,3)),(c_23,c_ij(2,3)),(s_13,s_ij(1,3)),(s_23,s_ij(2,3)),(d_cp,delta_cp)])
            print("Calculation for "+a+" to "+b+" started!")
            for i,l in enumerate(dist):
                eq_with_L=Equation.subs([(D_m_12,7.5*10**(-5)),(D_m_23,2.427*10**(-3)),(D_m_13,2.427*10**(-3)+7.5*10**(-5)),(L,l)])
                for j,en in enumerate(eneg_range):
                    Total_prob[i][j]=eq_with_L.subs([(E_nu,en)])
                del eq_with_L
            
            if a==e and b==mu:
                title=r"$P(\bar{\nu}_e\rightarrow\bar{\nu}_\mu)$"
            elif a==e and b==e:
                title=r"$P(\bar{\nu}_e\rightarrow\bar{\nu}_e)$"
            elif a==mu and b==mu:
                title=r"$P(\bar{\nu}_\mu\rightarrow\bar{\nu}_\mu)$"
            
            list_of_all.append([title,Total_prob])
            del Total_prob,Equation,title

        if bool and a!=b:
            bool=False

print("Probability calcualtions are done!")

print("Plots are starting...")
name_help=["e_to_e","e_to_mu","mu_to_mu"]
for events,name in zip(list_of_all,name_help):
    #events[0] is the title
    #events[1] is the Probability matrix

    plt.figure(figsize=(8,6))
    sns.heatmap(events[1],annot=False)

    #we are setting the axis
    ymin,ymax=plt.gca().get_ylim()
    xmin,xmax=plt.gca().get_xlim()  

    custom_ticks_y=np.linspace(ymin,ymax,9,dtype=int)
    custom_ticklabels_y=np.linspace(13000,0,9,dtype=int)
    custom_ticks_x=np.linspace(xmin,xmax,4)
    custom_ticklabels_x=[0.1,1,10,100]

    plt.gca().set_xticks(custom_ticks_x)
    plt.gca().set_xticklabels(custom_ticklabels_x)
    plt.gca().set_yticks(custom_ticks_y)
    plt.gca().set_yticklabels(custom_ticklabels_y)
    plt.xlabel("E(GeV)",fontsize=14)
    plt.ylabel("L(km)",fontsize=14)

    plt.xlim(right=int(n*0.89))

    plt.title(events[0],fontsize=14)
      
    plt.savefig('heatmaps/Hm_anti_Prob_%s_to.pdf' % (name))