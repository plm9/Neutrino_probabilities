import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from tqdm import tqdm
import seaborn as sns

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


def PMNS_param_matrix():
    m=Matrix([
        [c_12*c_13,s_12*c_13,s_13*exp(-d_cp*1j)],
        [-s_12*c_23-c_12*s_23*s_13*exp(d_cp*1j),c_12*c_23-s_12*s_23*s_13*exp(d_cp*1j),s_23*c_13],
        [s_12*s_23-c_12*c_23*s_13*exp(d_cp*1j),-c_12*s_23-s_12*c_23*s_13*exp(d_cp*1j),c_23*c_13]])

    return m

U=PMNS_param_matrix()

def flavor_to_index(a,b):
    if a == "e":
        index_a=0
    elif a == "mu":
        index_a=1
    elif a == "tau":
        index_a=2
    else:
        raise ValueError("Non existing neutrino flavor for initial state.")
    
    if b == "e":
        index_b=0
    elif b == "mu":
        index_b=1
    elif b == "tau":
        index_b=2
    else:
        raise ValueError("Non existing neutrino flavor for final state.")
    
    return index_a,index_b

def Prob_a_to_b(a,b):
    """ 
    Parameters:
    a: flavor of the neutrino in the initial state
    b: flavor of the neutrino in the final state

    Function:
    Calculates the probability of a neutrino of flavor a to transform into flavor b, as
    a function of the mixing angles ,the mass differences, the distance traveled and the 
    energy of the neutrino.
    """
    index_a,index_b=flavor_to_index(a,b)


    f=U[index_b,0]*conjugate(U[index_a,0])*U[index_a,1]*conjugate(U[index_b,1])
    s=U[index_b,0]*conjugate(U[index_a,0])*U[index_a,2]*conjugate(U[index_b,2])
    t=U[index_b,1]*conjugate(U[index_a,1])*U[index_a,2]*conjugate(U[index_b,2])

    first=re(f)*sin((D_m_21*L)/(4*E_nu))**2
    second=re(s)*sin((D_m_31*L)/(4*E_nu))**2
    third=re(t)*sin((D_m_32*L)/(4*E_nu))**2

    if a==b:
        return 1-4*(first+second+third)
    else:
        return -4*(first+second+third)

def Prob_anti_a_to_anti_b(a,b):
    index_a,index_b=flavor_to_index(a,b)

    f=conjugate(U[index_b,0])*U[index_a,0]*conjugate(U[index_a,1])*U[index_b,1]
    s=conjugate(U[index_b,0])*U[index_a,0]*conjugate(U[index_a,2])*U[index_b,2]
    t=conjugate(U[index_b,1])*U[index_a,1]*conjugate(U[index_a,2])*U[index_b,2]

    first=re(f)*sin((D_m_21*L)/(4*E_nu))**2
    second=re(s)*sin((D_m_31*L)/(4*E_nu))**2
    third=re(t)*sin((D_m_32*L)/(4*E_nu))**2

    if a==b:
        return 1-4*(first+second+third)
    else:
        return -4*(first+second+third)

def D_mass_param(i,j):
    if i==j:
        return ValueError("Probably 0 but not sure yet...")
    elif i==1:
        if j==2:
            return -D_m_21
        elif j==3:
            return -D_m_31
    elif i==2 :
        if j==1:
            return D_m_21
        elif j==3:
            return -D_m_32
    elif i==3:
        if j==1:
            return D_m_31
        elif j==2:            
            return D_m_32

def norm(z):
    return sqrt(re(z)**2+im(z)**2)

def Prob_a_to_b_alt(a,b):
    index_a,index_b=flavor_to_index(a,b)

    fst_rhs=0
    U=PMNS_param_matrix()
    for i in [0,1,2]:
        fst_rhs+=(norm(U[index_b,i])**2)* norm(U[index_a,i])**2

    scd_rhs=0
    for i in [0,1,2]:
        for j in [0,1,2]:
            if i>j:
                scd_rhs+=U[index_b,i]*conjugate(U[index_a,i])*conjugate(U[index_b,j])*U[index_a,j]*cos((D_mass_param(i+1,j+1))*L/(2*E_nu))

    return fst_rhs+(2*re(scd_rhs))
    


""" 
These functions are returning the sin and the cos respectively of the mixing angles.
"""
#The values where taken from the JUNO physics book
def s_ij(i,j):
    if (i==1 and j==2) or (i==2 and j==1):
        return 0.302**(1/2)
    elif (i==1 and j==3) or (i==3 and j==1):
        return 0.0227**(1/2)
    elif (i==2 and j==3) or (i==3 and j==2):
        return 0.413**(1/2)
    elif (i==j):
        raise ValueError("Not mixing! i=j")
    else:
        raise ValueError("Not excisting eigenstates! Must be 1,2 or 3.")

def c_ij(i,j):
    return (1-(s_ij(i,j)**2))**(1/2)

delta_cp=np.pi/2

#Defining the possible flavors as strings
e="e"
mu="mu"
tau="tau"

flavors=[e,mu]

#These variables are used fot the plots
n=10
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
            Equation=Prob_a_to_b_alt(a,b).subs([(s_12,s_ij(1,2)),(c_12,c_ij(1,2)),(c_13,c_ij(1,3)),(c_23,c_ij(2,3)),(s_13,s_ij(1,3)),(s_23,s_ij(2,3)),(d_cp,delta_cp)])
            #Equation_anti=Prob_anti_a_to_anti_b(a,b).subs([(s_12,s_ij(1,2)),(c_12,c_ij(1,2)),(c_13,c_ij(1,3)),(c_23,c_ij(2,3)),(s_13,s_ij(1,3)),(s_23,s_ij(2,3)),(d_cp,delta_cp)])
            print("Calculation for "+a+" to "+b+" started!")
            for i,l in enumerate(dist):
                eq_with_L=Equation.subs([(D_m_21,7.5*10**(-5)),(D_m_32,2.427*10**(-3)),(D_m_31,2.427*10**(-3)+7.5*10**(-5)),(L,l)])
                #eq_with_L_anti=Equation_anti.subs([(D_m_21,7.5*10**(-5)),(D_m_32,2.427*10**(-3)),(D_m_31,2.427*10**(-3)+7.5*10**(-5)),(L,l)])
                for j,en in enumerate(eneg_range):
                    Total_prob[i][j]=eq_with_L.subs([(E_nu,en)])#-eq_with_L_anti.subs([(E_nu,en)])
                del eq_with_L
            
            if a==e and b==mu:
                title=r"$P(\nu_e\rightarrow\nu_\mu)$"
            elif a==e and b==e:
                title=r"$P(\nu_e\rightarrow\nu_e)$"
            elif a==mu and b==mu:
                title=r"$P(\nu_\mu\rightarrow\nu_\mu)$"
            
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

    plt.title(events[0],fontsize=14) 
    plt.show()
    #plt.savefig('heatmaps/Hm_Prob_%s_to.pdf' % (name))
    
    


    