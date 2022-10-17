from turtle import title
import numpy as np
import package_func.my_functions as mf
import matplotlib.pyplot as plt
import package_func.Prob_funcs as pf
from sympy import rotations

#Constants
R_earth=6371 #km

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

def Long(cosZenith):
    return abs(-2*R_earth*cosZenith)

def CP_difference_with_35(a,b,cz,En):
    index_a,index_b=mf.flavor_to_index(a,b)
    term=0
    U=PMNS_matrix()
    for i in [0,1,2]:
        for j in [0,1,2]:
            if i>j:
                term+=(U[index_b,i]*U[index_a,i].conjugate()*U[index_b,j].conjugate()*U[index_a,j]).imag * np.sin(1.27*mf.D_mass(i+1,j+1)*Long(cz)/En)
    
    return 4*term 

def CP_difference_with_P_cp(cz,En):
    
    J=np.sin(2*theta_12)*np.sin(2*theta_23)*np.sin(2*theta_13)*np.cos(theta_13)*np.sin(d_cp)
    P_cp=J*np.sin(1.27*mf.D_mass(2,1)*Long(cz)/En)*np.sin(1.27*mf.D_mass(3,1)*Long(cz)/En)**2
    return P_cp


def main():
    N=200
    Cz=np.linspace(-1,0,N+1)
    En=np.logspace(-1,2,N+1)

    flavors=["e","mu","tau"]

    neutrinos=[r"$\nu_e$",r"$\nu_\mu$",r"$\nu_\tau$"]
    anti_neutrinos=[r"$\bar{\nu}_e$",r"$\bar{\nu}_\mu$",r"$\bar{\nu}_\tau$"]
    nu_index=0
    for a in flavors:
        data_list=[]
        titles=[]

        for b in flavors:
            #Prob=np.empty([N,N])
            #Prob_cp=np.empty([N,N])

            #Total_prob=np.empty([N,N])    
            #Total_prob_anti=np.empty([N,N])

            MSW_prob=np.empty([N,N])    
            MSW_prob_anti=np.empty([N,N])

            for i in range(N):
                for j in range(N):
                
                    Cs=0.5*(Cz[i]+Cz[i+1])
                    Eng=0.5*(En[j]+En[j+1])
                    
                    MSW_prob[i,j]=pf.Prob_a_to_b_Gen_MSW(a,b,Cs,Eng,ordering="NO")    
                    MSW_prob_anti[i,j]=pf.Prob_a_to_b_Gen_MSW(a,b,Cs,Eng,type="anti",ordering="NO")

                    #Total_prob[i,j]=pf.Prob_calc_General(a,b,Cs,Eng)
                    #Total_prob_anti[i,j]=pf.Prob_calc_General(a,b,Cs,Eng,type="anti")

                    #Prob[i,j]=CP_difference_with_35("e","mu",Cs,Eng)
                    #Prob_cp[i,j]=CP_difference_with_P_cp(Cs,Eng)
            data_list.append(MSW_prob)
            data_list.append(MSW_prob_anti)
            titles.append(a+"->"+b)


    #data_list=[MSW_prob_anti,MSW_prob,MSW_prob_anti-MSW_prob,Total_prob_anti,Total_prob,Total_prob_anti-Total_prob]#Prob[MSW_prob_anti,Prob_cp,Total_prob,MSW_prob]#Prob
    #titles=["MSW_anti","MSW","MSW anti-normal","General_anti","General","General anti-normal"]#"Eq_35""P_cp"
    

        fig, ax = plt.subplots(2, 3,figsize=(16, 18),sharex="col",sharey="row")
        indexer=0
        title_index=0
        for i in range(3):
            for j in range(2):
                data = data_list[indexer]
                im = ax[j, i].pcolormesh(En,Cz,data,cmap="hot",vmin=0,vmax=1)
                #plt.colorbar(im, ax=ax[j, i])
                ax[j,i].set_xscale("log")
                if indexer%2==0:
                    ax[j,i].set_title(titles[title_index])
                    title_index+=1  
                indexer+=1
        #plt.colorbar(im, ax=ax[j, i])
        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
        fig.colorbar(im, cax=cbar_ax) 
        ax[0,0].set_ylabel(neutrinos[nu_index],fontsize=16)
        ax[1,0].set_ylabel(anti_neutrinos[nu_index],fontsize=16)
        fig.text(0.07,0.5,r"$\cos(\theta)$",va='center',rotation="vertical",fontsize=16)
        fig.text(0.93,0.5,r"$P(\nu)$",fontsize=16)#va='center',rotation="vertical",
        fig.text(0.5, 0.04, r"$E_\nu$(GeV)", ha='center',fontsize=16)

        plt.savefig("Oscillation_probs_{}.pdf".format(a)) #I WILL CHANGE A AND THEN I LL ADD MSW ON THE NAME
        #plt.show()

        ## CP violation plots
        fig1, ax1 = plt.subplots(1, 3,figsize=(10.5, 4.5),sharex="col",sharey="row")
        indexer_CP=0
        title_index=0
        #for i in range(1):
        for j in range(3):
            data = data_list[indexer_CP]-data_list[indexer_CP+1]
            im = ax1[j].pcolormesh(En,Cz,data,cmap="bwr",vmin=-0.5,vmax=0.5)
            #plt.colorbar(im, ax=ax1[j])
            ax1[j].set_xscale("log")
            if indexer_CP%2==0:
                ax1[j].set_title(titles[title_index])
                title_index+=1  
            indexer_CP+=2
        fig1.subplots_adjust(right=0.85)
        cbar_ax = fig1.add_axes([0.87, 0.15, 0.02, 0.7])
        fig1.colorbar(im, cax=cbar_ax) 
        ax1[0].set_ylabel(neutrinos[nu_index]+"-"+anti_neutrinos[nu_index],fontsize=14)
        fig1.text(0.01,0.5,r"$\cos(\theta)$",va='center',rotation="vertical",fontsize=14)
        fig1.text(0.93,0.5,r"$P(\nu)-P(\bar{\nu})$",va='center',rotation="vertical",fontsize=14)
        fig1.text(0.5, 0.02, r"$E_\nu$(GeV)", ha='center',fontsize=14)

        plt.savefig("neutrinos-anti_{}.pdf".format(a))#I WILL CHANGE A AND THEN I LL ADD MSW ON THE NAME
        #plt.show()
        nu_index+=1
        del data_list,titles
        #break


if __name__=="__main__":
    main()
