import numpy as np
import matplotlib.pyplot as plt
from sympy import *
import seaborn as sns
#from Prob_calc import PMNS_param_matrix,flavor_to_index
from my_functions import *
from Prob_functions import *

#Defining symbols
s_12=Symbol("s_12",real=True)
s_23=Symbol("s_23",real=True)
s_13=Symbol("s_13",real=True)

c_12=Symbol("c_12",real=True)
c_23=Symbol("c_23",real=True)
c_13=Symbol("c_13",real=True)

d_cp=Symbol("d_cp",real=True)

D_m_21=Symbol("D_m_12",real=True)
D_m_31=Symbol("D_m_13",real=True)
D_m_32=Symbol("D_m_23",real=True)

L=Symbol("L",real=True,positive=True)
E_nu=Symbol("E_nu",real=True)

e="e"
mu="mu"
tau="tau"

flavors=[e,mu]

n=300
dist=np.linspace(0,13000,n)
eneg_range=np.logspace(-1,2,n)
Total_prob=np.zeros([n,n])

saving_list=[]

bool=True
for a in flavors:
    for b in flavors:
        if bool or a==b:
            expression=Prob_a_to_b_alt(a,b).subs([(s_12,sin(theta(1,2))),(c_12,cos(theta(1,2))),(c_13,cos(theta(1,3))),(c_23,cos(theta(2,3))),(s_13,sin(theta(1,3))),(s_23,sin(theta(2,3))),(d_cp,delta_cp),(D_m_21,D_mass(2,1)),(D_m_32,D_mass(3,2)),(D_m_31,D_mass(3,1))])
            Total_prob=np.zeros([n,n])
            print("Calculation for "+a+" to "+b+" started!")
            for i,l in enumerate(dist):
                eq_with_L=expression.subs([(L,l)])
                for j,en in enumerate(eneg_range):
                    Total_prob[i][j]=eq_with_L.subs([(E_nu,en)])
                del eq_with_L

            if a==e and b==mu:
                title=r"$P(\nu_e\rightarrow\nu_\mu)$"
            elif a==e and b==e:
                title=r"$P(\nu_e\rightarrow\nu_e)$"
            elif a==mu and b==mu:
                title=r"$P(\nu_\mu\rightarrow\nu_\mu)$"

            saving_list.append([title,Total_prob])
            del Total_prob,title,expression

        if bool and a!=b:
            bool=False

name_help=["e_to_e","e_to_mu","mu_to_mu"]
print("Plots are starting . . .")
for events,name in zip(saving_list,name_help):
    #events[0] is the title
    #events[1] is the Probability matrix

    plt.figure(figsize=(8,6))
    sns.heatmap(events[1],annot=False,cbar=False)
    sns.heatmap(events[1]).figure.axes[-1].set_ylabel(events[0],size=14) 

    #we are setting the axis
    ymin,ymax=plt.gca().get_ylim()
    xmin,xmax=plt.gca().get_xlim()  

    custom_ticks_y=np.linspace(ymin,ymax,6,dtype=int)
    custom_ticklabels_y=np.round(np.linspace(-1,0,6,dtype=float),2)
    custom_ticks_x=np.linspace(xmin,xmax,4)
    custom_ticklabels_x=[0.1,1,10,100]

    plt.gca().set_xticks(custom_ticks_x)
    plt.gca().set_xticklabels(custom_ticklabels_x)
    plt.gca().set_yticks(custom_ticks_y)
    plt.gca().set_yticklabels(custom_ticklabels_y)
    plt.xlabel("E(GeV)",fontsize=14)
    plt.ylabel(r"cos($\theta$)",fontsize=14)

    plt.xlim(right=int(n*0.89))

    #plt.title(events[0],fontsize=14) 
    plt.savefig('heatmaps/Me_prob_%s_to.pdf' % (name))
    plt.show()
    

print("Plots done!!!")

from Hardcode_plots import *

print("Difference Plots are starting . . .")
for events,hc_event,name in zip(saving_list,[prob_ne_ne,prob_ne_nmu,prob_nu_nmu],name_help):
    sns.heatmap(events[1]-hc_event,cbar=False)
    sns.heatmap(events[1]-hc_event).figure.axes[-1].set_ylabel("Difference of approaches",size=14) 

    #we are setting the axis
    ymin,ymax=plt.gca().get_ylim()
    xmin,xmax=plt.gca().get_xlim()  

    custom_ticks_y=np.linspace(ymin,ymax,6,dtype=int)
    custom_ticklabels_y=np.round(np.linspace(-1,0,6,dtype=float),2)
    custom_ticks_x=np.linspace(xmin,xmax,4)
    custom_ticklabels_x=[0.1,1,10,100]

    plt.gca().set_xticks(custom_ticks_x)
    plt.gca().set_xticklabels(custom_ticklabels_x)
    plt.gca().set_yticks(custom_ticks_y)
    plt.gca().set_yticklabels(custom_ticklabels_y)
    plt.xlabel("E(GeV)",fontsize=14)
    plt.ylabel(r"cos($\theta$)",fontsize=14)
    plt.title(events[0],fontsize=14)
    
    plt.xlim(right=int(n*0.89))
    plt.savefig("heatmaps/Difference_%s.pdf" %(name))
    plt.show()

print("Everything done")