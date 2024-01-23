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
    '2520000':'tttW_selected_Events.root',
    '2530000':'07E85166-4C72-3B4E-9572-2D070F2EF385.root',
    '260000':'07350478-4397-054C-B6F8-9F6979FF7357.root'
})
datasets.setdefault('tttJ', {}).update({
    'func': lambda dset, ptype: dset[ptype + 'tttJ'],
    'nrepos':3, 'repos':'2530000,260000,270000',
    '2530000':'46BC5096-1D58-1A45-81F3-F3E18CD4695A.root',
    '260000':'4668BDC5-735F-D14C-8B2B-5DC437CA5B25.root',
    '270000':'11CCD6D3-B6C1-8C45-9B20-08BF279D6BA4.root',
})
datasets.setdefault('tttt', {}).update({
    'func': lambda dset, ptype: dset[ptype + 'tttt'],
    'nrepos':1, 'repos':'2520000',
    '2520000':'4t_selected_Events.root',
})
tttW_style = {'histtype': 'step', 'density': False, 'lw': 1, 'zorder': 2,'label': 'tttW'}
for process,prcs in datasets.items():
    cosphi_0s=[]
    cosphi=[]
    for i in range(prcs['nrepos']):
        repo=(prcs['repos']).split(',')[i]
        Nfiles=len((prcs[repo]).split(','))
        for f in range(Nfiles):
            file=(prcs[repo].split(','))[f] #For now we only have 1 file in each repo but ok
            filepath='selected_Events/'+str(process)+'/'+str(repo)+'/'+str(file)
            cosphi_with_0s,cosphi_without_0s=Get_cosphi(filepath,process)
            cosphi_0s=np.concatenate([cosphi_0s,cosphi_with_0s])
            cosphi=np.concatenate([cosphi,cosphi_without_0s])
    print("Get_cosphi successfull for " + str(process))
    #Plot individual process distributions
    fig, ax = plt.subplots()
    hep.cms.label("", ax=ax, loc=0)
    ax.hist(cosphi,8,**tttW_style)
    ax.set_xlabel('Gen. level $\cos{\phi}$')
    ax.set_ylabel('Number of Events')
    fig.savefig('plots/'+str(process)+'_cosphi.pdf', bbox_inches='tight')
    