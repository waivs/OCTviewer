#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 13:51:18 2020

@author: jdrogers
"""


# %% imports
import numpy as np
import zipfile as zip
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt # we can drop this once we are done testing, this will only read the file and display will be handled by the calling script
from pathlib import Path
from scipy.ndimage import gaussian_filter

def readOCT(filename=None):
    '''
    read in Thorlabs OCT file and return intensity data
    '''
    if filename is None: 
        from tkinter.filedialog import askopenfilename # import inside if so if we never need this, it isn't loaded. If we are going to resort to dialog, speed of loading this matters not
        filename = askopenfilename()
    
    file = Path(filename)
    
    # check that file is a zipfile and Header.xml is present
    # TODO: error handling if the file we select is not a zip file, also check if it is a Thorlabs OCT file
    if zip.is_zipfile(file):
        print('we\'re good: selected file is a zipfile')
    else:
        print('uh oh, this looks like its not a zipfile and therefore not a Thorlabs OCT file')
    #with zip.ZipFile(file,'r') as myzip:
    #    print(myzip.infolist())    
    
    # open archive, list files, read Header.xml, determine size of array, read the image data
    # Note: this may break if the OCT mode is not 3D
    with zip.ZipFile(file,'r') as myzip:
        # files=myzip.namelist() # list the files contained in the archive
        # print(files)
        with myzip.open('Header.xml') as header:
            octinfo = header.read()
            root = ET.fromstring(octinfo)
        dims = root.findall("./Image/SizePixel/") # find the dimensions
        # dims is a list of dimensions with number of pixels in each 
        for dim in dims:
            print(f'{dim.tag}: {dim.text}')
            
        volsize = root.findall("./Image/SizeReal/")
            
        with myzip.open(r'data\Intensity.data') as intensitydata:
            data = (np.frombuffer(intensitydata.read(),dtype=np.float32).reshape(int(dims[2].text),int(dims[1].text),int(dims[0].text))).T
        
        # get xml data
        soctinfo = octinfo.decode() # convert bytes to string
        note = soctinfo[soctinfo.find('<Comment>')+9:soctinfo.find('</Comment>')]

    return(data, volsize, note)

    

# %% loop through all files
filepath = Path(r"./")
filematch = r'*_Mode3D.oct'
filenames = sorted(filepath.glob(filematch)) # sort by alphabet gives confocal, 2, 3, 4, 5

for file in filenames:
    # % Pick file and read in data
    #file = 'thisIsNotaZip.oct' # for testing a non-zip file OCT file
    # file = '/home/jdrogers/Documents/Data/OCT/Sajdak-BioScaffolds_0041_Mode3D.oct'
    imdata, volsize, note = readOCT(file)
    
    volsizemm = np.zeros(np.shape(volsize)[0]) # si
    for dim in np.arange(volsizemm.shape[0]):
        volsizemm[dim] = float(volsize[dim].text)
    
# %%            
    # % show image just for kicks
    # fig, ax = plt.subplots(nrows=2,ncols=2,num=0,figsize=(10,10),dpi=150, constrained_layout=True); fig.clear() # create, clear, create again ensure we always start fresh
    # fig, ax = plt.subplots(nrows=2,ncols=2,num=0);
    # fig.tight_layout()

    # fig, ax = plt.subplot_mosaic([['xy', 'yz'],['xz', '3d']],num = 0,figsize=(5, 5))
    # fig, ax = plt.subplot_mosaic([['xy', 'yz'],['xz', '3d']],
    #                           num = 0,
    #                           figsize=(5, 5),
    #                           width_ratios=(2, 1), height_ratios=(2, 1),
    #                           layout='constrained')
    
    
    # Create a Figure, which doesn't have to be square.
    fig = plt.figure(layout='constrained', num=0);fig.clear()
    fig = plt.figure(layout='constrained', num=0)
    ax0 = fig.add_subplot()
    ax0.set_aspect('equal')

    ax1 = ax0.inset_axes([1.05, 0, 0.5, 1], sharey=ax0)
    ax2 = ax0.inset_axes([0, -.55, 1, 0.5], sharex=ax0)
    ax3 = ax0.inset_axes([1.05,-.55,.5,.5], projection='3d')
    
    
    zstart = 150 #125 # start data at this plane to avoid artifacts
    zend = zstart+256 #512
    # thresh = 20 # the 'surface' is the first pixel along z to exceed this value
    thresh = imdata[zstart:zend].mean()+2*imdata[zstart:zend].std() # this should be more robust by choosing the threshold  that is some number of std's above the mean
    surfmap = (np.argmax(imdata[zstart:zend,:,:]>thresh,axis=0)) # surface map
    
    # given the surface map, set the colorbar limits to prevent a few low or high pixels from shifting the range too much
    nstd = 3 # number of std's away from mean to show
    hmin = surfmap.mean()-nstd*surfmap.std()
    hmax = surfmap.mean()+nstd*surfmap.std()
    surfmap = (np.argmax(imdata[zstart+int(hmin):zstart+int(hmax),:,:]>thresh,axis=0)) # surface map
    hmin = surfmap.min(); hmax=surfmap.max()
    
    ax0.imshow(imdata[zstart:zend].max(axis=0), cmap='gray',extent=[0,volsizemm[1],0,volsizemm[2]]); 
    ax0.title.set_text('maximum intensity projection');
    ax0.tick_params(labelbottom=False)
    ax0.set_ylabel('mm')
    ax1.imshow(imdata[zstart:zend,:,256].T, cmap='gray',extent=[zstart*(volsizemm[0]/imdata.shape[0]),zend*(volsizemm[0]/imdata.shape[0]),0,volsizemm[2]]); 
    ax1.title.set_text('y-z plane'); 
    ax1.tick_params(labelleft=False)
    ax1.annotate(note,(.8,.85),xycoords='subfigure fraction',color='white')
    ax2.imshow(imdata[zstart:zend,256,:], cmap='gray',extent=[0,volsizemm[1],zstart*(volsizemm[0]/imdata.shape[0]),zend*(volsizemm[0]/imdata.shape[0])]); 
    ax2.title.set_text('x-z plane'); 
    ax2.set_xlabel('mm');ax2.set_ylabel('mm')
    # ax[1,1].imshow(surfmap,cmap='viridis_r',vmin=hmin,vmax=hmax,extent=[0,volsizemm[1],0,volsizemm[2]]); ax[1,1].title.set_text('Surface'); ax[1,1].tick_params(labelleft=False)
    # ax[1,1].set_xlabel('mm')
    # ax[1,1].plot(imdata.mean(axis=(1,2)))
    fig.canvas.draw()
    
    #fig.savefig('out/'+file.name[:-4]+'.png')
    
    # %%
    if 1:
        # 3D perspective of surface
        fig1, ax1 = plt.subplots(subplot_kw={"projection": "3d"}, num=1,figsize=(10,10),dpi=150, constrained_layout=True); fig1.clear(); 
        fig1, ax1 = plt.subplots(subplot_kw={"projection": "3d"}, num=1,figsize=(10,10),dpi=150, constrained_layout=True);
        X = np.arange(imdata.shape[2])*volsizemm[2]/imdata.shape[2]
        Y = np.arange(imdata.shape[1])*volsizemm[1]/imdata.shape[1]
        Y,X=np.meshgrid(Y,X)
        Z = gaussian_filter(-surfmap*volsizemm[0]/imdata.shape[0],1.0)
        stride = 10
        # surf = ax1.plot_surface(X,Y,Z,cmap='viridis',rstride=stride,cstride=stride)
        surf = ax1.plot_surface(X,Y,Z,cmap='viridis',linewidth=0,edgecolor=None,antialiased='False',rstride=stride,cstride=stride)
        # surf = ax1.plot_surface(X,Y,Z,linewidth=0,antialiased='False',rstride=stride,cstride=stride)
        ax1.set_box_aspect(aspect=(1,1,0.1))
        ax1.set_aspect('equal') # only works in 3d for matplotlib > 3.6
        ax1.azim=10;ax1.elev=30
        ax1.set_zticks([Z.min(),Z.max()])
        ax1.set_zticklabels([f'{Z.min():.3f}',f'{Z.max():.3f}'])

        surf = ax3.plot_surface(X,Y,Z,cmap='viridis',linewidth=0,edgecolor=None,antialiased='False',rstride=stride,cstride=stride)
        # surf = ax1.plot_surface(X,Y,Z,linewidth=0,antialiased='False',rstride=stride,cstride=stride)
        ax3.set_box_aspect(aspect=(1,1,0.1))
        ax3.set_aspect('equal') # only works in 3d for matplotlib > 3.6
        ax3.azim=10;ax3.elev=30
        ax3.set_zticks([])
        ax3.set_xticks([])
        # ax3.set_zticklabels(None)


        fig1.canvas.draw()
        fig.canvas.draw()
        # fig1.savefig('out/'+file.name[:-4]+'3d'+'.png',bbox_inches='tight')
    
