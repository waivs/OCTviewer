# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 21:16:37 2022

@author: J.D. Rogers jeremy.rogers@wisc.edu
"""

# %% import
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from oct_converter.readers import e2e

# %% load data
#file = "/PATH/TO?FILE/filename.e2e"

try: # if file is defined above, use that, otherwise ask to choose the file with dialog 
    file
except NameError: # Note: tkinter askopen seems to be broken on mac
    from tkinter.filedialog import askopenfile, askdirectory # import inside so if you never need this, it isn't loaded. If we are going to resort to dialog, speed of loading doesn't matter
    # filepath = Path(askdirectory())
    f = askopenfile()
    file = f.name
    f.close()

# filematch = "*.E2E"
# filenames = sorted(filepath.glob(filematch)) 
# file = filenames[0].as_posix()

octdata = e2e.E2E(file)

# example from oct_converters library
#oct_volume.peek() # plots a montage of the volume
#oct_volume.save('fds_testing.avi')  # save volume as a movie
#oct_volume.save('fds_testing.png')  # save volume as a set of sequential images, fds_testing_[1...N].png
#oct_volume.save_projection('projection.png') # save 2D projection


# %% grab fundus images and OCT volumes
fundus_images = octdata.read_fundus_image()  # returns a Fundus image with additional metadata if available
fundus_ims = []
for imswmeta in fundus_images:
    fundus_ims.append(imswmeta.image)
print(f'Found {len(fundus_images)} fundus images')

oct_volumes = octdata.read_oct_volume()  # returns an OCT volume with additional metadata if available
oct_vols = []
for volwmeta in oct_volumes:
    oct_volume = volwmeta.volume
    octcube = np.concatenate(oct_volume,axis=0).reshape(-1,oct_volume[0].shape[0],oct_volume[0].shape[1])
    octcube = np.roll(octcube,1,axis=0) # account for bug in the file format or the reader where the first slice is actually stored at the end
    oct_vols.append(octcube)
    
print(f'Found {len(oct_volumes)} oct volumes')

# %% define an OCT viewer
def remove_keymap_conflicts(new_keys_set):
    for prop in plt.rcParams:
        if prop.startswith('keymap.'):
            keys = plt.rcParams[prop]
            remove_list = set(keys) & new_keys_set
            for key in remove_list:
                keys.remove(key)

# logscale = False # start with default
def oct_viewer(octvols,fundusims, num=None):
    
    if num==None:
        fig = plt.figure(figsize=(10,10),dpi=150,facecolor='gray')
        fig.clear() 
        if hasattr(fig, 'cid'): fig.canvas.mpl_disconnect(fig.cid) # if the canvas is already connected, disconnect to prevent multiple connections that lead to skipping slices on single keypress
        fig = plt.figure(num=fig.number,figsize=(10,10),dpi=150,facecolor='gray',layout='compressed')
        gs = GridSpec(2, 3,left=0.01, right=0.99, wspace=0.01, top=0.99, bottom=0.01, hspace=0.01)
        fig.axes[0] = fig.add_subplot(gs[0, 0:2]);fig.axes[0].axis('off')
        fig.axes[1] = fig.add_subplot(gs[1, :]);fig.axes[1].axis('off')
        fig.axes[2] = fig.add_subplot(gs[0, 2]);fig.axes[2].axis('off')
    else:
        fig = plt.figure(num=num,figsize=(10,10),dpi=150,facecolor='gray')
        fig.clear(); 
        if hasattr(fig, 'cid'): fig.canvas.mpl_disconnect(fig.cid) # if the canvas is already connected, disconnect to prevent multiple connections that lead to skipping slices on single keypress
        fig = plt.figure(num=num,figsize=(10,10),dpi=150,facecolor='gray',layout='compressed')
        gs = GridSpec(2, 3,left=0.01, right=0.99, wspace=0.01, top=0.99, bottom=0.01, hspace=0.01)
        fig.axes[0] = fig.add_subplot(gs[0, 0:2]);fig.axes[0].axis('off')
        fig.axes[1] = fig.add_subplot(gs[1, :]);fig.axes[1].axis('off')
        fig.axes[2] = fig.add_subplot(gs[0, 2]);fig.axes[2].axis('off')

    fig.vols = octvols
    fig.ims = fundusims
    # print(ax[1].volume)
    fig.volindex = 0
    fig.sliceindex = int(fig.vols[fig.volindex].shape[0]/2)
    fig.t = fig.axes[2].text(.01,.5,
                     f'next/prev slice: up/down \nnext/prev volume: left/right \n\nvolume: {fig.volindex+1} / {len(fig.vols)} \nslice: {fig.sliceindex}'
                     ,va="center", ha="left")

    fig.axes[1].imshow(fig.vols[fig.volindex][fig.sliceindex],cmap='gray')
    # fig.axes[1].set_title(f'volume {fig.volindex}; slice {fig.sliceindex}')
    fig.axes[0].imshow(fig.ims[fig.volindex],cmap='gray')
    fig.line = fig.axes[0].plot(np.array([256, 1536-256]),1536-256-fig.sliceindex*np.array([1024/97, 1024/97]),'g')
    fig.cid = fig.canvas.mpl_connect('key_press_event', process_key)
    return fig

def process_key(event):
    fig = event.canvas.figure
    if event.key == 'down': #'j':
        previous_slice(fig)
    elif event.key == 'up': #'k':
        next_slice(fig)
    elif event.key == 'left': #'k':
        previous_vol(fig)
    elif event.key == 'right': #'k':
        next_vol(fig)
        
    fig.t.set_text(f'next/prev slice: up/down \nnext/prev volume: left/right \n\nvolume: {fig.volindex+1} / {len(fig.vols)} \nslice: {fig.sliceindex}')
    fig.canvas.draw()
    
def previous_slice(fig):
    fig.sliceindex = (fig.sliceindex - 1) % fig.vols[fig.volindex].shape[0]
    # fig.axes[1].set_title(f'vol {fig.volindex}, slice {fig.sliceindex}')
    fig.axes[1].images[0].set_array(fig.vols[fig.volindex][fig.sliceindex])
    fig.axes[0].lines[-1].remove(); fig.line=fig.axes[0].plot(np.array([256, 1536-256]),1536-256-fig.sliceindex*np.array([1024/97, 1024/97]),'g')

def next_slice(fig):
    fig.sliceindex = (fig.sliceindex + 1) % fig.vols[fig.volindex].shape[0]
    # fig.axes[1].set_title(f'vol {fig.volindex}, slice {fig.sliceindex}')
    fig.axes[1].images[0].set_array(fig.vols[fig.volindex][fig.sliceindex])
    fig.axes[0].lines[-1].remove(); fig.line=fig.axes[0].plot(np.array([256, 1536-256]),1536-256-fig.sliceindex*np.array([1024/97, 1024/97]),'g')

def previous_vol(fig):
    if len(fig.vols)>1: fig.volindex = (fig.volindex - 1) % len(fig.vols)
    # fig.axes[1].set_title(f'vol {fig.volindex}, slice {fig.sliceindex}')
    fig.axes[1].images[0].set_array(fig.vols[fig.volindex][fig.sliceindex])
    fig.axes[0].images[0].set_array(fig.ims[fig.volindex])

def next_vol(fig):
    if len(fig.vols)>1: fig.volindex = (fig.volindex + 1) % len(fig.vols)
    # fig.axes[1].set_title(f'vol {fig.volindex}, slice {fig.sliceindex}')
    fig.axes[1].images[0].set_array(fig.vols[fig.volindex][fig.sliceindex])
    fig.axes[0].images[0].set_array(fig.ims[fig.volindex])



# %%
f0 = oct_viewer(oct_vols,fundus_ims,num=0)


 # %% run the viewer to display OCT slices and fundus photo


 # funduslateralities =[] # should be able to get the laterality of each image from either fundus or OCT volume metadata, but an apparent bug returns the wrong values for OCT, so use fundus for no
 # for ii in np.arange(len(fundus_images)): 
 #     funduslateralities.append(fundus_images[ii].laterality)
 # octlateralities =[] # should be able to get the laterality of each image from either fundus or OCT volume metadata, but an apparent bug returns the wrong values for OCT, so use fundus for no
 # for ii in np.arange(len(oct_volumes)): 
 #     octlateralities.append(oct_volumes[ii].laterality)

 # funduslaterality_index = 0 #funduslateralities.index(eye)
 # octlaterality_index = 0 #octlateralities.index(eye)
 # fundus_image = fundus_images[funduslaterality_index].image
 # oct_volume = oct_volumes[octlaterality_index].volume
 # octcube = np.concatenate(oct_volume,axis=0).reshape(-1,oct_volume[0].shape[0],oct_volume[0].shape[1])


 #a=oct_volume.contours.get('contour0'); b=np.concatenate(a,axis=0).reshape(len(a),-1); fig=plt.figure(num=1); fig.clear(); plt.imshow(b,aspect=b.shape[1]/b.shape[0]);fig.canvas.draw()
                             
 # # %% display in a figure
 # if 0:
 #     fig, ax = plt.subplots(nrows=2,ncols=1,num=0,figsize=(10,10),dpi=150); fig.clear();
 #     fig, ax = plt.subplots(nrows=2,ncols=1,num=0,figsize=(10,10),dpi=150)
 #     plt.subplots_adjust(left=0.01, bottom=0.01, right=0.98, top=0.97, wspace=0.01, hspace=0.06)

 #     # Could do this in a loop, but this makes it clear which PMTs get which labels adn scaling
 #     ax[0].imshow(fundus_image.image, cmap='gray')
 #     ax[0].set_title('')
 #     ax[0].axis('off')

 #     linenum = 47
 #     ax[1].imshow((oct_volume.volume[linenum]),cmap='gray', vmax=0.6*oct_volume.volume[linenum].max())
 #     ax[1].set_title('')
 #     ax[1].axis('off')

 #     fig.canvas.draw()
