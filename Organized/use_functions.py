import numpy as np
from my_functions import *

#Constants
R_earth=6,371 #km

#angles
theta_12=33.45 #deg
theta_23=42.1
theta_13=8.62

d_cp=90 #deg -> this can change as much as we want

s_12=np.sin(theta_12)
s_23=np.sin(theta_23)
s_13=np.sin(theta_13)

c_12=np.cos(theta_12)
c_23=np.cos(theta_23)
c_13=np.cos(theta_13)


def PMNS_matrix():
    m=np.array([
        [c_12*c_13,s_12*c_13,s_13*np.exp(-d_cp*1j)],
        [-s_12*c_23-c_12*s_23*s_13*np.exp(d_cp*1j),c_12*c_23-s_12*s_23*s_13*np.exp(d_cp*1j),s_23*c_13],
        [s_12*s_23-c_12*c_23*s_13*np.exp(d_cp*1j),-c_12*s_23-s_12*c_23*s_13*np.exp(d_cp*1j),c_23*c_13]])

    return m

def L(cosZenith):
    return -2*R_earth*cosZenith

def CP_difference_with_35(a,b,cz,En):
    index_a,index_b=flavor_to_index(a,b)
    term=0
    U=PMNS_matrix()
    for i in [0,1,2]:
        for j in [0,1,2]:
            term+=(U[index_b,i]*U[index_a,i].conjugate()*U[index_b,j].conjugate()*U[index_a,j]).imag * np.sin(1.27*D_mass(i+1,j+1)*cz/En)
    
    return 4*term 

def main():
    CC=-1
    En=10
    print(CP_difference_with_35("e","mu",CC,En))

if __name__=="__main__":
    main()
