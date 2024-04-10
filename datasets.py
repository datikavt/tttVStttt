import numpy as np

def get_datasets():
    datasets={} #initialize an empty dictionary to store the datasets that we want to use for plotting
    datasets.setdefault('tttW', {}).update({
        'func': lambda dset, ptype: dset[ptype + 'tttW'],
        'nrepos':5, 'repos':'2520000,2530000,260000,231218_110037,local', #repos separated by commas
        # 'nrepos':3, 'repos':'2520000,2530000,260000',
        '2520000':'9660E5EF-E305-7546-B26E-9A759DBC1976.root',
        '2530000':'07E85166-4C72-3B4E-9572-2D070F2EF385.root',
        '260000':'07350478-4397-054C-B6F8-9F6979FF7357.root,ECACC29E-698C-EE42-A460-BAD444727AFA.root', #Filenames should be separated by commas
        '231218_110037':'NanoAOD_1.root,NanoAOD_2.root,NanoAOD_3.root,NanoAOD_4.root',
        'local':'NanoAOD_1.root',
    })
    # More tttJ datasets: https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset+primary_dataset%3DTTTJ*
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
        '2520000':'1323716F-65CC-1B46-96FE-2E660766A235.root,2D36D54E-8CCB-E144-90C5-3312F0805C7F.root,3876D6A8-8D88-424D-8DE2-408C95A1AF78.root,5F6252A3-1634-404F-A493-2A36345FDA1B.root,6AAF523E-FCEE-ED42-B161-C6204614BDF8.root,9135C82E-46A9-0B4A-AE88-2149604E80F8.root,AED740DE-A0C3-9248-846F-3B6EA9B689EE.root',
        # '2520000':'1323716F-65CC-1B46-96FE-2E660766A235.root,2D36D54E-8CCB-E144-90C5-3312F0805C7F.root,5F6252A3-1634-404F-A493-2A36345FDA1B.root'
    })
    return datasets