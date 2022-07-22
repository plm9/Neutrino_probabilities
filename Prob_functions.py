import numpy as np
from sympy import *

import package_func.my_functions as mf

L=Symbol("L",real=True,positive=True)
E_nu=Symbol("E_nu",real=True)

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
    index_a,index_b=mf.flavor_to_index(a,b)

    U=mf.PMNS_param_matrix()

    f=U[index_b,0]*conjugate(U[index_a,0])*U[index_a,1]*conjugate(U[index_b,1])
    s=U[index_b,0]*conjugate(U[index_a,0])*U[index_a,2]*conjugate(U[index_b,2])
    t=U[index_b,1]*conjugate(U[index_a,1])*U[index_a,2]*conjugate(U[index_b,2])

    first=re(f)*sin((mf.D_m_21*L)/(4*E_nu))**2
    second=re(s)*sin((mf.D_m_31*L)/(4*E_nu))**2
    third=re(t)*sin((mf.D_m_32*L)/(4*E_nu))**2

    if a==b:
        return 1-4*(first+second+third)
    else:
        return -4*(first+second+third)

def Prob_anti_a_to_anti_b(a,b):
    index_a,index_b=mf.flavor_to_index(a,b)

    U=mf.PMNS_param_matrix()

    f=conjugate(U[index_b,0])*U[index_a,0]*conjugate(U[index_a,1])*U[index_b,1]
    s=conjugate(U[index_b,0])*U[index_a,0]*conjugate(U[index_a,2])*U[index_b,2]
    t=conjugate(U[index_b,1])*U[index_a,1]*conjugate(U[index_a,2])*U[index_b,2]

    first=re(f)*sin((mf.D_m_21*L)/(4*E_nu))**2
    second=re(s)*sin((mf.D_m_31*L)/(4*E_nu))**2
    third=re(t)*sin((mf.D_m_32*L)/(4*E_nu))**2

    if a==b:
        return 1-4*(first+second+third)
    else:
        return -4*(first+second+third)

def Prob_a_to_b_alt(a,b):
    index_a,index_b=mf.flavor_to_index(a,b)
    
    U=mf.PMNS_param_matrix()

    fst_rhs=0
    U=mf.PMNS_param_matrix()
    for i in [0,1,2]:
        fst_rhs+=(mf.norm(U[index_b,i])**2)* mf.norm(U[index_a,i])**2

    scd_rhs=0
    for i in [0,1,2]:
        for j in [0,1,2]:
            if i>j:
                scd_rhs+=U[index_b,i]*conjugate(U[index_a,i])*conjugate(U[index_b,j])*U[index_a,j]*cos((mf.D_mass_param(i+1,j+1))*L/(2*E_nu))

    return fst_rhs+(2*re(scd_rhs))


def Prob_a_to_b_General(a,b,type):
    index_a,index_b=mf.flavor_to_index(a,b)

    if type == "anti":
        U=conjugate(mf.PMNS_param_matrix())
    else:
        U=mf.PMNS_param_matrix()

    re_sum=0
    im_sum=0
    for i in [0,1,2]:
        for j in [0,1,2]:
            if i>j:
                term=U[index_b,i]*conjugate(U[index_a,i])*conjugate(U[index_b,j])*U[index_a,j]
                re_sum+=re(term)*(sin((mf.D_mass_param(i+1,j+1))*L/(4*E_nu)))**2
                im_sum+=im(term)*(sin((mf.D_mass_param(i+1,j+1))*L/(2*E_nu)))
                del term
    
    return mf.delta_Kro(a,b)-4*re_sum-2*im_sum

def Prob_a_to_b_Gen_MSW(a,b,type="normal"):
    index_a,index_b=mf.flavor_to_index(a,b)

    if type == "anti":
        U=conjugate(mf.PMNS_param_matrix())
    else:
        U=mf.PMNS_param_matrix()

    re_sum=0
    im_sum=0
    for i in [0,1,2]:
        for j in [0,1,2]:
            if i>j:
                if i+1==1 and j+1==3:
                    mass_term=mf.MSW_Dmass(i+1,j+1)
                else:
                    mass_term=mf.D_mass_param(i+1,j+1)
                term=U[index_b,i]*conjugate(U[index_a,i])*conjugate(U[index_b,j])*U[index_a,j]
                re_sum+=re(term)*(sin((mass_term)*L/(4*E_nu)))**2
                im_sum+=im(term)*(sin((mass_term)*L/(2*E_nu)))
                del term
    returner=mf.delta_Kro(a,b)-4*re_sum-2*im_sum
    return returner

