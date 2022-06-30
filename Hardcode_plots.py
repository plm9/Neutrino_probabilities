from my_functions import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

delta_cp=np.pi/2

J=(1/8)*np.sin(2*deg_to_rad(theta_12))*np.sin(2*deg_to_rad(theta_23))*np.sin(2*deg_to_rad(theta_13))*np.cos(deg_to_rad(theta_13))*np.sin(delta_cp)

P_cp=8*J*sin((D_mass(2,1)*L)/(4*E_nu))*sin((D_mass(3,1)*L)/(4*E_nu))**2

S_23=sin((D_mass(3,2)*L)/(4*E_nu))**2
S_12=sin((D_mass(2,1)*L)/(4*E_nu))**2

P_ne_nmu=(np.sin(theta_23)**2)*(np.sin(2*theta_13)**2)*S_23+(np.cos(theta_23)**2)*(np.sin(2*theta_12)**2)*S_12-P_cp
P_ne_ne=1-(np.sin(2*theta_13)**2)*S_23-(np.cos(theta_13)**4)*(np.sin(2*theta_12)**2)*S_12
P_nmu_nmu=1-(4*(np.cos(theta_13)**2)*(np.sin(theta_23)**2)*(1-(np.cos(theta_13)**2)*(np.sin(theta_23)**2))*S_23)-((np.cos(theta_23)**2)*(np.sin(2*theta_12)**2)*S_12)

n=300
dist=np.linspace(0,13000,n)
eneg_range=np.logspace(-1,2,n)

prob_ne_ne=np.zeros([n,n])
prob_ne_nmu=np.zeros([n,n])
prob_nu_nmu=np.zeros([n,n])


for i,l in enumerate(dist):
    eq_ee=P_ne_ne.subs([(L,l)])
    eq_em=P_ne_nmu.subs([(L,l)])
    eq_mm=P_nmu_nmu.subs([(L,l)])
    for j,en in enumerate(eneg_range):
        prob_ne_ne[i][j]=eq_ee.subs([(E_nu,en)])
        prob_ne_nmu[i][j]=eq_em.subs([(E_nu,en)])
        prob_nu_nmu[i][j]=eq_mm.subs([(E_nu,en)])
    del eq_ee,eq_mm,eq_em


plt.figure(figsize=(8,6))
sns.heatmap(prob_ne_ne,cbar=False)
sns.heatmap(prob_ne_ne).figure.axes[-1].set_ylabel(r"$P(\nu_{e}\rightarrow\nu_{e})$",size=14) 
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
plt.ylabel(r"$cos(\theta)$",fontsize=14)

plt.xlim(right=int(n*0.89))
#plt.title(r"$P(\nu_{e}\rightarrow\nu_{e})$",fontsize=14) 
plt.savefig("heatmaps/heatmaps_from_paper/HM_prob_e_e.pdf")
plt.show()

plt.figure(figsize=(8,6))
sns.heatmap(prob_ne_nmu,cbar=False)
sns.heatmap(prob_ne_nmu).figure.axes[-1].set_ylabel(r"$P(\nu_{e}\rightarrow\nu_{\mu})$",size=14) 
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
plt.ylabel(r"$cos(\theta)$",fontsize=14)

plt.xlim(right=int(n*0.89))
#plt.title(r"$P(\nu_{e}\rightarrow\nu_{\mu})$",fontsize=14) 
plt.savefig("heatmaps/heatmaps_from_paper/HM_prob_e_mu.pdf")
plt.show()

plt.figure(figsize=(8,6))
sns.heatmap(prob_nu_nmu,cbar=False)
sns.heatmap(prob_nu_nmu).figure.axes[-1].set_ylabel(r"$P(\nu_{\mu}\rightarrow\nu_{\mu})$",size=14) 
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
plt.ylabel(r"$cos(\theta)$",fontsize=14)

plt.xlim(right=int(n*0.89))
#plt.title(r"$P(\nu_{\mu}\rightarrow\nu_{\mu})$",fontsize=14) 
plt.savefig("heatmaps/heatmaps_from_paper/HM_prob_mu_mu.pdf")
plt.show()
