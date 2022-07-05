from my_functions import *
import matplotlib.pyplot as plt
import seaborn as sns

from Prob_functions import *

e="e"
mu="mu"
tau="tau"

flavors=[e,mu]

n=200
dist=np.linspace(0,13000,n)
eneg_range=np.logspace(-1,2,n)

saving_list=[]
saving_list_anti=[]

saving_list_diff=[]

bool=True
for a in flavors:
    for b in flavors:
        if bool or a==b:
            expression=Prob_a_to_b_General(a,b,"normal").subs([(s_12,sin(theta(1,2))),(c_12,cos(theta(1,2))),(c_13,cos(theta(1,3))),(c_23,cos(theta(2,3))),(s_13,sin(theta(1,3))),(s_23,sin(theta(2,3))),(d_cp,delta_cp),(D_m_21,D_mass(2,1)),(D_m_32,D_mass(3,2)),(D_m_31,D_mass(3,1))])
            expression_anti=Prob_a_to_b_General(a,b,"anti").subs([(s_12,sin(theta(1,2))),(c_12,cos(theta(1,2))),(c_13,cos(theta(1,3))),(c_23,cos(theta(2,3))),(s_13,sin(theta(1,3))),(s_23,sin(theta(2,3))),(d_cp,delta_cp),(D_m_21,D_mass(2,1)),(D_m_32,D_mass(3,2)),(D_m_31,D_mass(3,1))])
            expression_diff=(Prob_a_to_b_General(a,b,"normal")-Prob_a_to_b_General(a,b,"anti")).subs([(s_12,sin(theta(1,2))),(c_12,cos(theta(1,2))),(c_13,cos(theta(1,3))),(c_23,cos(theta(2,3))),(s_13,sin(theta(1,3))),(s_23,sin(theta(2,3))),(d_cp,delta_cp),(D_m_21,D_mass(2,1)),(D_m_32,D_mass(3,2)),(D_m_31,D_mass(3,1))])
            Total_prob=np.zeros([n,n])
            Total_prob_anti=np.zeros([n,n])
            Total_prob_diff=np.zeros([n,n])
            print("Calculation for "+a+" to "+b+" started!")
            for i,l in enumerate(dist):
                eq_with_L=expression.subs([(L,l)])
                eq_with_L_a=expression_anti.subs([(L,l)])
                eq_with_L_diff=expression_diff.subs([(L,l)])
                for j,en in enumerate(eneg_range):
                    Total_prob[i][j]=eq_with_L.subs([(E_nu,en)])
                    Total_prob_anti[i][j]=eq_with_L_a.subs([(E_nu,en)])
                    Total_prob_diff[i][j]=eq_with_L_diff.subs([(E_nu,en)])
                del eq_with_L,eq_with_L_a,eq_with_L_diff

            if a==e and b==mu:
                title=r"$P(\nu_e\rightarrow\nu_\mu)$"
                title_anti=r"$P(\bar{\nu}_e\rightarrow\bar{\nu}_\mu)$"
            elif a==e and b==e:
                title=r"$P(\nu_e\rightarrow\nu_e)$"
                title_anti=r"$P(\bar{\nu}_e\rightarrow\bar{\nu}_e)$"
            elif a==mu and b==mu:
                title=r"$P(\nu_\mu\rightarrow\nu_\mu)$"
                title_anti=r"$P(\bar{\nu}_\mu\rightarrow\bar{\nu}_\mu)$"

            title_diff="CP-violation"

            saving_list.append([title,Total_prob])
            saving_list_anti.append([title_anti,Total_prob_anti])
            saving_list_diff.append([title_diff,Total_prob_diff])
            del Total_prob,title,expression

        if bool and a!=b:
            bool=False

name_help=["e_to_e","e_to_mu","mu_to_mu"]
print("Plots are starting . . .")
for events,antievents,diff_events,name in zip(saving_list,saving_list_anti,saving_list_diff,name_help):
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
    #plt.savefig('heatmaps/general/Gen_prob_%s_to.pdf' % (name))
    plt.show()

    #Here we plot for the anti neutinos
    plt.figure(figsize=(8,6))
    sns.heatmap(antievents[1],annot=False,cbar=False)
    sns.heatmap(antievents[1]).figure.axes[-1].set_ylabel(antievents[0],size=14) 

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
    #plt.savefig('heatmaps/general/Me_prob_anti_%s_to.pdf' % (name))
    plt.show()


    plt.figure(figsize=(8,6))
    sns.heatmap(diff_events[1],cbar=False,annot=False)
    sns.heatmap(diff_events[1]).figure.axes[-1].set_ylabel(events[0]+"-"+antievents[0],size=14) 
    
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

    plt.title(diff_events[0],fontsize=14) 
    #plt.savefig('heatmaps/general/Gen_diff_3pi2_%s_to.pdf' % (name))
    plt.show()


print("Plots done!!!")