from ast import Expression
import numpy as np
from sympy import *
from Prob_gen import Total_prob

from my_functions import  *

def Calculator(flavors,Prob_func_expr,nb):
    dist=np.linspace(0,13000,nb)
    eneg_range=np.logspace(-1,2,nb)

    saving_list=[]

    for a in flavors:
        for b in flavors:
                expression=Prob_func_expr.subs([(s_12,sin(theta(1,2))),(c_12,cos(theta(1,2))),(c_13,cos(theta(1,3))),(c_23,cos(theta(2,3))),(s_13,sin(theta(1,3))),(s_23,sin(theta(2,3))),(d_cp,delta_cp),(D_m_21,D_mass(2,1)),(D_m_32,D_mass(3,2)),(D_m_31,D_mass(3,1))])
                Total_prob=np.zeros([nb,nb])
                for i,l in enumerate(dist):
                    eq_with_L=expression.subs([(L,l)])
                    for j,en in enumerate(eneg_range):
                        Total_prob[i][j]=eq_with_L.subs([(E_nu,en)])
                
                    del eq_with_L
                
                if a=="e" and b=="mu":
                    title=r"$P(\nu_e\rightarrow\nu_\mu)$"
                    title_anti=r"$P(\bar{\nu}_e\rightarrow\bar{\nu}_\mu)$"
                elif a=="e" and b=="tau":
                    title=r"$P(\nu_e\rightarrow\nu_\tau)$"
                    title_anti=r"$P(\bar{\nu}_e\rightarrow\bar{\nu}_\tau)$"
                elif a=="mu" and b=="tau":
                    title=r"$P(\nu_\mu\rightarrow\nu_\tau)$"
                    title_anti=r"$P(\bar{\nu}_\mu\rightarrow\bar{\nu}_\tau)$"
                elif a=="e" and b=="e":
                    title=r"$P(\nu_e\rightarrow\nu_e)$"
                    title_anti=r"$P(\bar{\nu}_e\rightarrow\bar{\nu}_e)$"
                elif a=="mu" and b=="mu":
                    title=r"$P(\nu_\mu\rightarrow\nu_\mu)$"
                    title_anti=r"$P(\bar{\nu}_\mu\rightarrow\bar{\nu}_\mu)$"
                elif a=="tau" and b=="tau":
                    title=r"$P(\nu_\tau\rightarrow\nu_\tau)$"
                    title_anti=r"$P(\bar{\nu}_\tau\rightarrow\bar{\nu}_\tau)$"

