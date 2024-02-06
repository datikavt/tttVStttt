import mplhep as hep #CAT
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
hep.style.use("CMS") 
import cmsstyle as CMS # For ROOT plotting
from Get_cosphi1 import Get_cosphi
from Get_cosphi1 import Get_delta_phi
from datasets import get_datasets
import argparse  # I/O
import math



parser = argparse.ArgumentParser(description='selecting events from a ROOT file.')
parser.add_argument('-obs','--observable',type=str, help='specify the observables that you want to plot')
parser.add_argument('-sel','--selection_type',type=bool,help='specify if you want to work with Mass_requirement selected or not',default=False)
args = parser.parse_args()
datasets=get_datasets() 
#Plot the D and cosphi distributions
if args.observable=='D':

    tttW_style = {'histtype': 'step', 'density': False, 'lw': 3, 'zorder': 2,'label': 'tttW','color':'#1fd216'}
    tttW_style_bar={'lw': 3, 'zorder': 2,'label': 'tttW'}
    tttJ_style={'histtype': 'step', 'density': False, 'lw': 3, 'zorder': 3,'label': 'tttJ','color':'orange'}
    tttJ_style_bar={'lw': 3, 'zorder': 3,'label': 'tttJ'}
    tttt_style={'histtype': 'step', 'density': True, 'lw': 3, 'zorder': 2,'label': 'tttt','color':'#159cb2'}
    ttt_style={'histtype':'step','density':True,'lw':2,'zorder':1,'label':'ttt','fill':False, 'color':'#bc5858'}

    cosphi_arrays = {}
    for process,prcs in datasets.items():
        cosphi=[]
        for i in range(prcs['nrepos']):
            repo=(prcs['repos']).split(',')[i]
            Nfiles=len((prcs[repo]).split(','))
            for f in range(Nfiles):
                file=(prcs[repo].split(','))[f] #For now we only have 1 file in each repo but ok
                if args.selection_type==True:
                    print(str(process) + "/"+str(repo))
                    filepath='selected_Events/'+str(process)+'/Mass_Pair_requirement/'+str(repo)+'/'+str(file)
                else:
                    filepath='selected_Events/'+str(process)+'/'+str(repo)+'/'+str(file)
                cosphi_file=Get_cosphi(filepath,process)
                cosphi=np.concatenate([cosphi,cosphi_file])
        print("Get_cosphi successfull for " + str(process))
        #Plot individual process distributions
        fig, (ax, axr) = plt.subplots(
            ncols=2,
            gridspec_kw=dict(width_ratios=[3, 0.25], wspace=0.5),
        )
        hep.cms.label("", ax=ax, loc=0)
        ax.set_xlabel('Gen. level $\cos{\phi}$')
        ax.set_ylabel('Number of Pairs')
        if process=='tttt':
            ax.set(ylim=(11000,15000),xlim=(-1,1))
            ax.hist(cosphi,8,histtype='step',density=False,lw=2,label='tttt',color='#159cb2')
        else:    
            ax.set(ylim=(600,1500),xlim=(-1,1))
            if process=='tttW':
                ax.hist(cosphi,8,**tttW_style)
            else:
                ax.hist(cosphi,8,**tttJ_style)
        ax.legend()
        D_proxy=-3*np.average(np.asarray(cosphi)) #This includes all the events even when tt pair was not found to have the same mother
        print("D for "+str(process)+" is " + str(D_proxy) )
        axr.set(ylim=(-0.4,0.1),ylabel="Gen. level D",xlim=(0,1))
        axr.scatter(0.5,D_proxy)
        axr.errorbar(0.5,D_proxy,xerr=0,yerr=(1/np.sqrt(len(cosphi))))
        axr.set_xticklabels([])
        if args.selection_type==True:
            fig.savefig('plots/Mass_Pair_requirement/'+str(process)+'_cosphi.pdf', bbox_inches='tight')
        else:
            fig.savefig('plots/'+str(process)+'_cosphi.pdf', bbox_inches='tight')
        cosphi_arrays[process]=cosphi



    ttt_cosphi=np.concatenate([cosphi_arrays['tttW'],cosphi_arrays['tttJ']])
    fig, (ax, axr) = plt.subplots(
          ncols=2,
            gridspec_kw=dict(width_ratios=[3, 0.5], wspace=0.5),
        )
    hep.cms.label("", ax=ax, loc=0)
    D_proxy_ttt=-3*np.average(np.asarray(ttt_cosphi))
    print("D for ttt is " + str(D_proxy_ttt))
    D_proxy_tttt=-3*np.average(np.asarray(cosphi_arrays['tttt']))
    ax.set_xlabel('Gen. level $\cos{\phi}$')
    ax.set_ylabel('Event density')
    ax.hist(cosphi_arrays['tttt'],8,**tttt_style)
    ax.hist(ttt_cosphi,8,**ttt_style)
    ax.set(ylim=(0.4,1),xlim=(-1,1))
    ax.legend()
    axr.set(ylim=(-0.4,0.1),ylabel="Gen. level D",xlim=(0,1))
    axr.scatter(0.2,D_proxy_ttt,marker='H',color='#bc5858')
    axr.errorbar(0.2,D_proxy_ttt,xerr=0,yerr=(1/np.sqrt(len(ttt_cosphi))),color='#bc5858')
    axr.scatter(0.4,D_proxy_tttt,color='#159cb2')
    axr.errorbar(0.4,D_proxy_tttt,xerr=0,yerr=(1/np.sqrt(len(cosphi_arrays['tttt']))),color='#159cb2')
    # axr.scatter(0.6,D_proxy_ttt,marker='H',color='#bc5858',label='D ttt')
    # axr.errorbar(0.6,D_proxy_ttt,xerr=0,yerr=(1/np.sqrt(len(ttt_cosphi))),color='#bc5858')
    # axr.scatter(0.8,D_proxy_tttt,color='#159cb2',label='D tttt')
    # axr.errorbar(0.8,D_proxy_tttt,xerr=0,yerr=(1/np.sqrt(len(cosphi_arrays['tttt']))),color='#159cb2')
    axr.set_xticklabels([]) # x doesn't need any ticks
    axr.legend()
    if args.selection_type==True:
        fig.savefig('plots/Mass_Pair_requirement/All_processes_cosphi.png', bbox_inches='tight')
    else:
        fig.savefig('plots/All_processes_cosphi.png', bbox_inches='tight')




#Get the A lab_cosphi and A delta phis
A_lab_cosphi_arrays = {}
for process,prcs in datasets.items():
    lab_cosphi=[]
    delta_phi=[] 
    for i in range(prcs['nrepos']):
        repo=(prcs['repos']).split(',')[i]
        Nfiles=len((prcs[repo]).split(','))
        for f in range(Nfiles):
            file=(prcs[repo].split(','))[f]
            if args.selection_type==True:
                filepath='selected_Events/'+str(process)+'/Mass_Pair_requirement/'+str(repo)+'/'+str(file)
            else:
                filepath='selected_Events/'+str(process)+'/'+str(repo)+'/'+str(file)
            lab_cosphi_file=Get_cosphi(filepath,process) #The A_lab_cos
            delta_phi_file=Get_delta_phi(filepath,process) #The A delta phi of the file
            lab_cosphi=np.concatenate([lab_cosphi,lab_cosphi_file])
            delta_phi=np.concatenate([delta_phi,delta_phi_file])
    print("Some of the delta phis are: " +str(delta_phi[:10]))
    cos_phi_lab_positive=np.sum(np.array(lab_cosphi) > 0)
    cos_phi_lab_negative=np.sum(np.array(lab_cosphi)<0)
    delta_phi_more=np.sum(np.array(delta_phi)>math.pi/2)
    delta_phi_less=np.sum(np.array(delta_phi)<math.pi/2)
    A_cos_phi_lab=(cos_phi_lab_positive-cos_phi_lab_negative)/(cos_phi_lab_positive+cos_phi_lab_negative)
    print(A_cos_phi_lab)
    A_delta_phi=(delta_phi_more-delta_phi_less)/(delta_phi_more+delta_phi_less)
    print("A_delta_phi is " + str(A_delta_phi))
    A_lab_cosphi_arrays[process+str('_A')]=np.array(A_cos_phi_lab)
    A_lab_cosphi_arrays[process+str('_A_delta_phi')]=np.array(A_delta_phi)
    A_lab_cosphi_arrays[process + str('_delta_phi')]=np.array(delta_phi)
    A_lab_cosphi_arrays[process + str('_cosphi')]=np.array(lab_cosphi) #Store the whole lab cosphi information not just A for later



ttt_lab_cosphi=np.concatenate([A_lab_cosphi_arrays['tttW_cosphi'],A_lab_cosphi_arrays['tttJ_cosphi']])
ttt_lab_cosphi_positive=np.sum(np.array(ttt_lab_cosphi) > 0)
ttt_lab_cosphi_negative = np.sum(np.array(ttt_lab_cosphi) <0)
A_lab_cosphi_ttt=(ttt_lab_cosphi_positive-ttt_lab_cosphi_negative)/(ttt_lab_cosphi_positive+ttt_lab_cosphi_negative)
print("A lab cosphi 3 top is " + str(A_lab_cosphi_ttt))
print("A lab cosphi 4 top is " + str(A_lab_cosphi_arrays['tttt_A']))


ttt_delta_phi=np.concatenate([A_lab_cosphi_arrays['tttW_delta_phi'],A_lab_cosphi_arrays['tttJ_delta_phi']])
ttt_delta_more=np.sum(np.array(ttt_delta_phi)>math.pi/2)
ttt_delta_less=np.sum(np.array(ttt_delta_phi)<math.pi/2)
A_delta_phi_ttt=(ttt_delta_more-ttt_delta_less)/(ttt_delta_more+ttt_delta_less)



#Plot the A's:
fig, axr=plt.subplots()
axr.set(ylim=(-0.25,0.25),ylabel="Gen. level Coefficients",xlim=(0,1))
x_points_ttt=np.array([0.2,0.6])
x_points_tttt=np.array([0.25,0.65])
x_points_ttbar=np.array([0.3,0.7])
y_points_ttt=np.array([A_lab_cosphi_ttt,A_delta_phi_ttt])
y_points_tttt=np.array([A_lab_cosphi_arrays['tttt_A'],A_lab_cosphi_arrays['tttt_A_delta_phi']])
y_points_ttbar=np.array([0.167080,0.103314]) #taken from https://www.hepdata.net/record/ins1742786



# colors=np.array(['#bc5858','#159cb2'])
yerrs=[(1/np.sqrt((len(ttt_lab_cosphi)))),(1/np.sqrt((len(A_lab_cosphi_arrays['tttt_cosphi'])))),(1/np.sqrt((len(ttt_delta_phi)))),1/np.sqrt((len(A_lab_cosphi_arrays['tttt_delta_phi'])))]
y_errs_ttbar=[0.003454,0.003201]

for i in range(len(x_points_ttt)):
    axr.scatter(x_points_ttt[i],y_points_ttt[i],marker='H',label='ttt',color='#bc5858')
    axr.errorbar(x_points_ttt[i],y_points_ttt[i],xerr=0,yerr=yerrs[0::2][i],ecolor='#bc5858')
    axr.scatter(x_points_tttt[i],y_points_tttt[i],marker='H',label='tttt',color='#159cb2')
    axr.errorbar(x_points_tttt[i],y_points_tttt[i],xerr=0,yerr=yerrs[1::2][i],ecolor='#159cb2')
    axr.scatter(x_points_ttbar[i],y_points_ttbar[i],marker='H',label='ttbar',color='#6a9253')
    axr.errorbar(x_points_ttbar[i],y_points_ttbar[i],xerr=0,yerr=y_errs_ttbar[i],ecolor='#6a9253')
# axr.scatter(0.2,A_lab_cosphi_ttt,marker='H',color='#bc5858',label='ttt $A^{lab}_{\cos \phi}$')
# axr.errorbar(0.2,A_lab_cosphi_ttt,xerr=0,yerr=(1/np.sqrt((len(ttt_lab_cosphi)))),color='#bc5858')
# axr.scatter(0.25,A_lab_cosphi_arrays['tttt_A'],marker='H',color='#159cb2',label='tttt $A^{lab}_{\cos \phi}$')
# axr.errorbar(0.25,A_lab_cosphi_arrays['tttt_A'],xerr=0,yerr=(1/np.sqrt((len(A_lab_cosphi_arrays['tttt_cosphi'])))),color='#159cb2')
# axr.scatter(0.6,A_delta_phi_ttt,marker='o',color='#bc5858',label='ttt $A_{|\Delta \phi_{ll}|}$')
# axr.errorbar(0.6,A_delta_phi_ttt,xerr=0,yerr=(1/np.sqrt((len(ttt_delta_phi)))),color='#bc5858')
# axr.scatter(0.65,A_lab_cosphi_arrays['tttt_A_delta_phi'],marker='o',color='#159cb2',label='tttt $A_{|\Delta \phi_{ll}|}$')
# axr.errorbar(0.65,A_lab_cosphi_arrays['tttt_A_delta_phi'], xerr=0,yerr=(1/np.sqrt((len(A_lab_cosphi_arrays['tttt_delta_phi'])))))
# axr.set_xticklabels([]) # x doesn't need any ticks
def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))

ticks=[0.225,0.625]
labels=["$A^{lab}_{\cos \phi}$","$A_{|\Delta \phi_{ll}|}$"]
axr.set_xticks(ticks=ticks)
axr.set_xticklabels(labels=labels)
legend_without_duplicate_labels(axr)
# axr.axhline(y=0, color='green', linestyle='--', label='Zero Line')
if args.selection_type==True:
    fig.savefig('plots/Mass_Pair_requirement/Coefficients.png', bbox_inches='tight')
else:
    fig.savefig('plots/Coefficients.png', bbox_inches='tight')





    