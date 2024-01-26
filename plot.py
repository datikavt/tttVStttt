import mplhep as hep #CAT
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
hep.style.use("CMS") 
import cmsstyle as CMS # For ROOT plotting
from Get_cosphi1 import Get_cosphi

datasets={} #initialize an empty dictionary to store the datasets that we want to use for plotting
datasets.setdefault('tttW', {}).update({
    'func': lambda dset, ptype: dset[ptype + 'tttW'],
    'nrepos':3, 'repos':'2520000,2530000,260000', #repos separated by commas
    '2520000':'9660E5EF-E305-7546-B26E-9A759DBC1976.root',
    '2530000':'07E85166-4C72-3B4E-9572-2D070F2EF385.root',
    '260000':'07350478-4397-054C-B6F8-9F6979FF7357.root,ECACC29E-698C-EE42-A460-BAD444727AFA.root',
})
datasets.setdefault('tttJ', {}).update({
    'func': lambda dset, ptype: dset[ptype + 'tttJ'],
    'nrepos':4, 'repos':'2520000,2530000,260000,270000',
    '2520000':'DBB31367-6197-CC4F-AE0C-CED640EE24A2.root',
    '2530000':'46BC5096-1D58-1A45-81F3-F3E18CD4695A.root',
    '260000':'4668BDC5-735F-D14C-8B2B-5DC437CA5B25.root',
    '270000':'11CCD6D3-B6C1-8C45-9B20-08BF279D6BA4.root',
})
datasets.setdefault('tttt', {}).update({
    'func': lambda dset, ptype: dset[ptype + 'tttt'],
    'nrepos':1, 'repos':'2520000',
    '2520000':'1323716F-65CC-1B46-96FE-2E660766A235.root,2D36D54E-8CCB-E144-90C5-3312F0805C7F.root,3876D6A8-8D88-424D-8DE2-408C95A1AF78.root,5F6252A3-1634-404F-A493-2A36345FDA1B.root,6AAF523E-FCEE-ED42-B161-C6204614BDF8.root,9135C82E-46A9-0B4A-AE88-2149604E80F8.root',
})

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
    axr.errorbar(0.5,D_proxy,xerr=0,yerr=(1/np.sqrt(len(cosphi)/2)))
    axr.set_xticklabels([])
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
# tttW_hist,bins=np.histogram(cosphi_arrays['tttW'],bins=8)
# tttJ_hist=np.histogram(cosphi_arrays['tttJ'],bins=8)[0]
# ax.hist(cosphi_arrays['tttW'],8,**tttW_style)
# ax.hist(cosphi_arrays['tttJ'],8,**tttJ_style)
ax.set_xlabel('Gen. level $\cos{\phi}$')
ax.set_ylabel('Event density')
ax.hist(cosphi_arrays['tttt'],8,**tttt_style)
# labels=['tttW','tttJ']
# colors=[tttW_style['color'],tttJ_style['color']] 
# ax.hist([cosphi_arrays['tttW'],cosphi_arrays['tttJ']],8,**ttt_style,label=labels,color=colors)
ax.hist(ttt_cosphi,8,**ttt_style)
ax.set(ylim=(0.4,1),xlim=(-1,1))
ax.legend()
axr.set(ylim=(-0.4,0.1),ylabel="Gen. level D",xlim=(0,1))
axr.scatter(0.25,D_proxy_ttt,marker='H',color='#bc5858')
axr.errorbar(0.25,D_proxy_ttt,xerr=0,yerr=(1/np.sqrt(len(ttt_cosphi)/2)),color='#bc5858')
axr.scatter(0.75,D_proxy_tttt,color='#159cb2')
axr.errorbar(0.75,D_proxy_tttt,xerr=0,yerr=(1/np.sqrt(len(cosphi_arrays['tttt'])/2)),color='#159cb2')
axr.set_xticklabels([]) # x doesn't need any ticks
fig.savefig('plots/All_processes_cosphi.pdf', bbox_inches='tight')

