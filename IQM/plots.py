## JAKE A TUTMAHER
## JOHNS HOPKINS UNIVERSITY
## MCQUEEN LABORATORY   
## THE INSTITUTE FOR QUANTUM MATTER
## DEPARTMENT OF PHYSICS, DEPARTMENT OF CHEMISTRY, DEPARTMENT OF MATERIALS SCIENCE AND ENGINEERING
##
## CONTACT: jtutmah1@jhu.edu
##################################################################################################

"""Plotting routines for VASP output files. They support spin-polarized and non-spin-polarized, as
   well as orbital character plots. Currently does not support Compounds with only p-electrons.
   Currently optimized for arrays generated by VASPread.
   
   ---
   
   eledosplot:           Provides a standard DOS plot - automatically recognizes spin-polarized input 
                         from VASPread.
   elebandplot:          Provides a standard bandplot - automatically recognizes spin-polarized input
                         from VASPread. Requires input of both the special points and the labels.
   elebanddosplot:       Provides a bandplot with an adjacent verticl DOS plot - automatically recog-
                         nizes spin-polarized input from VASPread. Requires input of the special pts
			 and the labels."""

####################################################################################################

#=================================  MODULES  =============================================
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from numpy import linalg as LA
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
#import IQM.vasp as vasp 

#===========================  ELECTRONIC DOS PLOT  =======================================
def eledosplot(dosx,dosy,fermilevel,dir_save):
    
    #Define Full Figure
    fig1 = plt.figure(figsize=(10,7),dpi=200)
    matplotlib.rc('xtick',labelsize=20)
    matplotlib.rc('ytick',labelsize=20)

    #Min and Max
    xmin=np.amin(dosx)
    xmax=np.amax(dosx)
    ymin=np.amin(dosy)
    ymax=np.amax(dosy)

    #Define DOS Plot Location and Attributes
    ax1 = plt.subplot()
    ax1.tick_params(axis='x',which='both',bottom='on', top='on',labelbottom='on')
    ax1.axes.get_yaxis().set_ticklabels([])
    ax1.axis([xmin,xmax,ymin,ymax])
    ax1.set_ylabel('Amplitude',size=20)
    ax1.set_xlabel('Energy (eV)',size=20)
    
    #CHECK FOR SPIN POLARIZED - PLOT DOS
    if dosy.shape[0]==2:
	   ax1.plot(dosx[0,:],dosy[0,:],color='darkred',linewidth=2)
	   ax1.plot(dosx[0,:],dosy[1,:],color='darkred',linewidth=2,linestyle=u'dashed')
    else:
	   ax1.plot(dosx[0,:],dosy[0,:],color='darkred',linewidth=2)
    
    #Plot Fermi Level
    ax1.vlines(fermilevel,ymin,ymax,colors='orange',linestyles=u'dashed',linewidth=2)
       
    plt.savefig(dir_save+'/ElectronicDOS.pdf')
    return

#=============================  SUMMARIZE BAND CHARACTERS  =======================================
def summarize_dos_characters(characters,total_dos):
    #DETERMINE DIMENSIONS OF CHARACTER ARRAY
    spin_length = characters.shape[0]
    val_length = characters.shape[1]
    char_length = characters.shape[2]
    
    #SET DIMENSIONS OF INTERNAL ARRAY
    spin_array = np.zeros((val_length,char_length))
    total_array = np.zeros((val_length))
    
    #COMBINE SPIN DATA
    if spin_length==2:
	  for i in range(spin_length):
	    total_array += total_dos[i,:]
	    spin_array += characters[i,:,:]
	#spin_array = np.sqrt(spin_array)
    else:
	   total_array = total_dos[0,:]
	   spin_array = characters[0,:,:]
    
    #SUMMARIZE CHARACTERS
    out_array = np.zeros((val_length,4))
    if char_length==9:
	  for i in range(val_length):
	    total = sum(spin_array[i,:])
	    if total==0:
		  out_array[i,0]=spin_array[i,0]
		  out_array[i,1]=(sum(spin_array[i,1:3]))
		  out_array[i,2]=(sum(spin_array[i,4:]))
	    else:
		  out_array[i,0]=spin_array[i,0]*(total_array[i]/total)
		  out_array[i,1]=(sum(spin_array[i,1:3]))*(total_array[i]/total)
		  out_array[i,2]=(sum(spin_array[i,4:]))*(total_array[i]/total)		
    elif char_length==16:
	  for i in range(val_length):
	    total = sum(spin_array[i,:])
	    if total==0:
		  out_array[i,0] = spin_array[i,0]
		  out_array[i,1] = (sum(spin_array[i,1:3]))
		  out_array[i,2] = (sum(spin_array[i,4:9]))
		  out_array[i,3] = (sum(spin_array[i,10:]))
	    else:
		  out_array[i,0] = spin_array[i,0]*(total_array[i]/total)
		  out_array[i,1] = (sum(spin_array[i,1:3]))*(total_array[i]/total)
		  out_array[i,2] = (sum(spin_array[i,4:9]))*(total_array[i]/total)
		  out_array[i,3] = (sum(spin_array[i,10:]))*(total_array[i]/total)		
    else: #LENGTH ALREADY EQUALS 4
	   out_array = spin_array
    
    return out_array
	    

#===========================  ELECTRONIC ORBITAL DOS PLOT  =======================================
def orbital_dosplot(dosx,dosy,total_dosy,fermilevel,dir_save):

    #Define Full Figure
    fig1 = plt.figure(figsize=(10,7),dpi=200)
    matplotlib.rc('xtick',labelsize=20)
    matplotlib.rc('ytick',labelsize=20)

    #Min and Max
    xmin=np.amin(dosx)
    xmax=np.amax(dosx)

    #Define DOS Plot Location and Attributes
    ax1 = plt.subplot()
    ax1.tick_params(axis='x',which='both',bottom='on', top='on',labelbottom='on')
    ax1.axes.get_yaxis().set_ticklabels([])
    ax1.set_ylabel('Amplitude',size=20)
    ax1.set_xlabel('Energy (eV)',size=20)

    #Plot DOS
    summarized_dosy = summarize_dos_characters(dosy,total_dosy)
    ymin=np.amin(summarized_dosy)
    ymax=np.amax(summarized_dosy)
    ax1.axis([xmin,xmax,ymin,ymax])
    
    ax1.plot(dosx[0,:],summarized_dosy[:,0],color='k',linewidth=2,label="s-dos")
    ax1.plot(dosx[0,:],summarized_dosy[:,1],color='r',linewidth=2,label="p-dos")
    ax1.plot(dosx[0,:],summarized_dosy[:,2],color='g',linewidth=2,label="d-dos")
    ax1.plot(dosx[0,:],summarized_dosy[:,3],color='b',linewidth=2,label="f-dos")

    #Plot Fermi Level
    ax1.vlines(fermilevel,ymin,ymax,colors='orange',linestyles=u'dashed',linewidth=2)
    plt.legend()

    plt.savefig(dir_save+'/Orbital_ElectronicDOS.pdf')
    return

#============================  ELECTRONIC BS PLOT  =======================================
def elebandplot(ymat,reclat,kpoints,fermilevel,kpt_spec,labels,dir_save):
	
	#GEN X-VALUES
    xmat = np.dot(reclat,kpoints)
	lengthx = len(xmat[0,:])
	xvalues = np.zeros(lengthx)
	total = 0
	for i in range(1,lengthx):
	    distance = LA.norm(xmat[:,i]-xmat[:,i-1])
	    total+=distance
	    xvalues[i] = total
	
	#Inputs (Minimums and Maximums)
	xmax = max(xvalues)
	ymin = fermilevel-10
	ymax = fermilevel+10
	#ymin = np.amin(ymat)
	#ymax = np.amax(ymat)
	
	#Define Full Figure
	width=10 #inch
	height=7 #inch
	fig1 = plt.figure(figsize=(width,height),dpi=200)
	matplotlib.rc('xtick',labelsize=20)
	matplotlib.rc('ytick',labelsize=20)
	
	#Define Band Plot Location and Attributions
	ax1 = plt.subplot()
	ax1.tick_params(axis='x',which='both',bottom='on', top='off',labelbottom='off')
	ax1.axis([0,xmax,ymin,ymax])
	ax1.set_ylabel('Energy (eV)',size=20)
	
	#Plot Electronic Bands
	if ymat.shape[0] == 2: #SPIN POLARIZED
	    for n in range(len(ymat[0,0,:])):
		  nup_plot, = ax1.plot(xvalues,ymat[0,:,n],color='k',linewidth=3,label='Spin Up')
		  ndown_plot, = ax1.plot(xvalues,ymat[1,:,n],color='gray',linestyle='--',linewidth=2,label='Spin Down')
	else:
	    for n in range(len(ymat[0,0,:])):
		  nplot, = ax1.plot(xvalues,ymat[0,:,n],color='k',linewidth=2,label='Energy')	
	
	#Set Dimensions of Image Canvas
	start=0.08
	finish=0.96
	scale = finish-start
	plt.subplots_adjust(left=start+0.01,right=finish,bottom=0.1,top=0.92)
	#Generate First and Last Label
	plt.figtext(start,0.05,labels[0],size=20)
	plt.figtext(start+scale,0.05,labels[len(labels)-1],size=20)
	#Set Remaining Partitions and Labels
	x_spec = np.dot(reclat,kpt_spec)
	distance=0
	for m in range(len(x_spec[0,:])-1):
	    interval = LA.norm(x_spec[:,m+1]-x_spec[:,m])
	    distance += interval
	    ax1.vlines(distance,ymin,ymax,colors='black',linewidth=2)
	    plt.figtext(start+scale*distance/max(xvalues),0.05,labels[m+1],size='20')

	#Plot Fermi Level
	fermi_plot = ax1.hlines(fermilevel,0,xmax,colors='orange',linestyles=u'dashed',linewidth=2,label='Fermi Level')
	
	#Legend
	if ymat.shape[0]==2:
	    plt.legend(handles=[nup_plot,ndown_plot,fermi_plot])
	else:
	    plt.legend(handles=[nplot,fermi_plot])	
	plt.savefig(dir_save+"/ElectronicBS.pdf")
	return

#=============================  SUMMARIZE BAND CHARACTERS  =======================================
def summarize_characters(characters):
    #DETERMINE DIMENSIONS OF CHARACTER ARRAY
    spin_length = characters.shape[0]
    kpt_length = characters.shape[1]
    band_length = characters.shape[2]
    char_length = characters.shape[3]
    
    #SET DIMENSIONS OF INTERNAL ARRAY
    spin_array = np.zeros((kpt_length,band_length,char_length))
    
    #COMBINE SPIN DATA
    if spin_length==2:
	   for i in range(spin_length):
	       spin_array += characters[i,:,:,:]**2
	       spin_array = np.sqrt(spin_array)
    else:
	   spin_array = characters[0,:,:,:]
    
    #SUMMARIZE CHARACTERS
    out_array = np.zeros((kpt_length,band_length,4))
    if char_length==9:
	for i in range(kpt_length):
	    for j in range(band_length):
		  tot = (spin_array[i,j,0]**2+sum(spin_array[i,j,1:3]**2)+sum(spin_array[i,j,4:]**2))
		  out_array[i,j,0]=spin_array[i,j,0]**2
		  out_array[i,j,1]=(sum(spin_array[i,j,1:3]**2))
		  out_array[i,j,2]=(sum(spin_array[i,j,4:]**2))
		  if tot !=0:
		    out_array[i,j,0]=out_array[i,j,0]/tot
		    out_array[i,j,1]=out_array[i,j,1]/tot
		    out_array[i,j,2]=out_array[i,j,2]/tot
		
    elif char_length==16:
	for i in range(kpt_length):
	    for j in range(band_length):
		  tot = (spin_array[i,j,0]**2+sum(spin_array[i,j,1:3]**2)+sum(spin_array[i,j,4:9]**2)+sum(spin_array[i,j,10:]**2))
		  out_array[i,j,0] = spin_array[i,j,0]**2
		  out_array[i,j,1] = (sum(spin_array[i,j,1:3]**2))
		  out_array[i,j,2] = (sum(spin_array[i,j,4:9]**2))
		  out_array[i,j,3] = (sum(spin_array[i,j,10:]**2))
		  if tot!=0:
		      out_array[i,j,0]=out_array[i,j,0]/tot
		      out_array[i,j,1]=out_array[i,j,1]/tot
		      out_array[i,j,2]=out_array[i,j,2]/tot
		      out_array[i,j,3]=out_array[i,j,3]/tot
		
    else: #LENGTH ALREADY EQUALS 4
	out_array = spin_array
    
    return out_array

#==========================  ELECTRONIC ORBITAL BAND PLOT ================================
def orbital_elebandplot(ymat,reclat,kpoints,fermilevel,kpt_spec,labels,characters,dir_save):
    """
    This method is partifularly usefule for orbital bandplots. See examples for more detail
    """
    
    #PLOT FOR EACH MULTICOLOR BAND
    def rgbline(ax, k, e, red, green, blue, alpha=1.):
	pts = np.array([k, e]).T.reshape(-1, 1, 2)
	seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
    
	nseg = len(k)-1
	r = [0.5*(red[p]+red[p+1]) for p in range(nseg)]
	g = [0.5*(green[p]+green[p+1]) for p in range(nseg)]
	b = [0.5*(blue[p]+blue[p+1]) for p in range(nseg)]
	a = np.ones(nseg, np.float)*alpha
	lc = LineCollection(seg, colors=zip(r,g,b,a), linewidth = 2)
	ax.add_collection(lc)

    #GEN X-VALUES
    xmat = np.dot(reclat,kpoints)
    lengthx = len(xmat[0,:])
    xvalues = np.zeros(lengthx)
    total = 0
    for i in range(1,lengthx):
	   distance = LA.norm(xmat[:,i]-xmat[:,i-1])
	   total+=distance
	   xvalues[i] = total

    #Inputs (Minimums and Maximums)
    xmax = max(xvalues)
    ymin = fermilevel-10
    ymax = fermilevel+10
    #ymin = np.amin(ymat)
    #ymax = np.amax(ymat)

    #Define Full Figure
    width=10 #inch
    height=7 #inch
    fig1 = plt.figure(figsize=(width,height),dpi=200)
    matplotlib.rc('xtick',labelsize=20)
    matplotlib.rc('ytick',labelsize=20)

    #Define Band Plot Location and Attributions
    ax1 = plt.subplot()
    ax1.tick_params(axis='x',which='both',bottom='on', top='off',labelbottom='off')
    ax1.axis([0,xmax,ymin,ymax])
    ax1.set_ylabel('Energy (eV)',size=20)
    
    #SUMMARIZE CHARACTERS
    characters_new = summarize_characters(characters)

    #Plot Electronic Bands
    if ymat.shape[0] == 2: #SPIN POLARIZED
	   for n in range(len(ymat[0,0,:])):
	       nup_plot = rgbline(ax1, xvalues, ymat[0,:,n], characters_new[:,n,1], characters_new[:,n,2], characters_new[:,n,3], alpha=1.)
	       ndown_plot = rgbline(ax1, xvalues, ymat[1,:,n], characters_new[:,n,1], characters_new[:,n,2], characters_new[:,n,3], alpha=1.)		    
    else:
	   for n in range(len(ymat[0,0,:])):	    
	       nplot = rgbline(ax1,xvalues,ymat[0,:,n],characters_new[:,n,1],characters_new[:,n,2],characters_new[:,n,3],alpha=1.)	

    #Set Dimensions of Image Canvas
    start=0.08
    finish=0.96
    scale = finish-start
    plt.subplots_adjust(left=start+0.01,right=finish,bottom=0.1,top=0.92)
    
    #Generate First and Last Label
    plt.figtext(start,0.05,labels[0],size=20)
    plt.figtext(start+scale,0.05,labels[len(labels)-1],size=20)
    
    #Set Remaining Partitions and Labels
    x_spec = np.dot(reclat,kpt_spec)
    distance=0
    for m in range(len(x_spec[0,:])-1):
	   interval = LA.norm(x_spec[:,m+1]-x_spec[:,m])
	   distance += interval
	   ax1.vlines(distance,ymin,ymax,colors='black',linewidth=2)
	   plt.figtext(start+scale*distance/max(xvalues),0.05,labels[m+1],size='20')

    #Plot Fermi Level
    fermi_plot = ax1.hlines(fermilevel,0,xmax,colors='orange',linestyles=u'dashed',linewidth=2,label='Fermi Level')
    
    #LEGEND
    splot = plt.hlines(-2000,0,xmax,color='k',label='s-band')
    pplot = plt.hlines(-2000,0,xmax,color='r',label='p-band')
    dplot = plt.hlines(-2000,0,xmax,color='g',label='d-band')
    fplot = plt.hlines(-2000,0,xmax,color='b',label='f-band')
    plt.legend(handles=[splot,pplot,dplot,fplot,fermi_plot])    
    
    plt.savefig(dir_save+"/Orbital_ElectronicBS.pdf")
    return
#==========================  ELECTRONIC BS AND DOS PLOT  =================================
def elebanddosplot(ymat,reclat,kpoints,fermilevel,kpt_spec,labels,dosx,dosy,dir_save, SpinPol=False):
	
    #GEN X-VALUES
    xmat = np.dot(reclat,kpoints)
    lengthx = len(xmat[0,:])
    xvalues = np.zeros(lengthx)
    total = 0
    for i in range(1,lengthx):
	   distance = LA.norm(xmat[:,i]-xmat[:,i-1])
	   total+=distance
	   xvalues[i] = total

    #Inputs (Minimums and Maximums)
    xmax = max(xvalues)
    ymin = fermilevel-10
    ymax = fermilevel+10
    #ymin = np.amin(ymat)
    #ymax = np.amax(ymat)

    #Define Full Figure
    width=14
    height=7.5
    fig1 = plt.figure(figsize=(width,height),dpi=200)
    gs = gridspec.GridSpec(2,7)
    matplotlib.rc('xtick',labelsize=20)
    matplotlib.rc('ytick',labelsize=20)

    #Define Band Plot Location and Attributions
    ax1 = plt.subplot(gs[:,0:5])
    ax1.tick_params(axis='x',which='both',bottom='on', top='off',labelbottom='off')
    ax1.axis([0,xmax,ymin,ymax])
    ax1.set_ylabel('Energy (eV)',size=20)

    #Plot Electronic Bands
    if ymat.shape[0]==2:
	   for n in range(len(ymat[0,0,:])):
	       nup_plot, = ax1.plot(xvalues,ymat[0,:,n],color='k',linewidth=3,label='Spin Up')
	       ndown_plot, = ax1.plot(xvalues,ymat[1,:,n],color='gray',linestyle='--',linewidth=2,label='Spin Down')
    else:
	   for n in range(len(ymat[0,0,:])):
	       nplot, = ax1.plot(xvalues,ymat[0,:,n],color='k',linewidth=2,label='Energy')    

    #Set Dimensions of Image Canvas
    start=0.06
    finish=0.65
    scale = finish-start
    plt.subplots_adjust(left=start+0.01,right=0.9,bottom=0.1,top=0.92)
    #Generate First and Last Label
    plt.figtext(start,0.05,labels[0],size=20)
    plt.figtext(start+scale,0.05,labels[len(labels)-1],size=20)
    #Set Remaining Partitions and Labels
    x_spec = np.dot(reclat,kpt_spec)
    distance=0
    for m in range(len(x_spec[0,:])-1):
	   interval = LA.norm(x_spec[:,m+1]-x_spec[:,m])
	   distance += interval
	   ax1.vlines(distance,ymin,ymax,colors='black',linewidth=2)
	   plt.figtext(start+scale*distance/max(xvalues),0.05,labels[m+1],size='20')

    #Plot Fermi Level
    fermi_plot = ax1.hlines(fermilevel,0,xmax,colors='orange',linestyles=u'dashed',linewidth=2,label='Fermi Level')

    #Define DOS Plot Location and Attributes
    ax2 = plt.subplot(gs[:,5:])
    ax2.tick_params(axis='x',which='both',bottom='on', top='on',labelbottom='off')
    ax2.axes.get_yaxis().set_ticklabels([])
    ax2.axis([0,1.2*np.amax(dosy),ymin,ymax])

    #Plot DOS
    if dosy.shape[0]==2:
	   dosplot_up, = ax2.plot(dosy[0,:],dosx[0,:],color='darkred',linewidth=2,label='DOS-Up')
	   dosplot_down, = ax2.plot(dosy[1,:],dosx[0,:],color='darkred',linewidth=2,linestyle='--',label='DOS-Down')
    else:
	   dosplot, = ax2.plot(dosy[0,:],dosx[0,:],color='darkred',linewidth=2,label='DOS')
    
    #Plot Fermi Level 
    ax2.hlines(fermilevel,0,1.2*np.amax(dosy),colors='orange',linestyles=u'dashed',linewidth=2)
    
    #Plot Legend
    #Legend
    if ymat.shape[0]==2:
	   plt.legend(handles=[nup_plot,ndown_plot,dosplot_up,dosplot_down,fermi_plot],bbox_to_anchor=(1.4,1.01)) 
    else:
	   plt.legend(handles=[nplot,dosplot,fermi_plot],bbox_to_anchor=(1.4,1.01))    
    
    #Plot Save
    plt.savefig(dir_save+"/ElectronicBSandDOS.pdf")
    return 

#==========================  ELECTRONIC BS AND DOS PLOT  =================================
def orbital_elebanddosplot(ymat,reclat,kpoints,fermilevel,kpt_spec,labels,characters,dosx,dosy,total_dosy,dir_save):
    
    #PLOT FOR EACH MULTICOLOR BAND
    def rgbline(ax, k, e, red, green, blue, alpha=1.):
	pts = np.array([k, e]).T.reshape(-1, 1, 2)
	seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

	nseg = len(k)-1
	r = [0.5*(red[p]+red[p+1]) for p in range(nseg)]
	g = [0.5*(green[p]+green[p+1]) for p in range(nseg)]
	b = [0.5*(blue[p]+blue[p+1]) for p in range(nseg)]
	a = np.ones(nseg, np.float)*alpha
	lc = LineCollection(seg, colors=zip(r,g,b,a), linewidth = 2)
	ax.add_collection(lc)
    
    #GEN X-VALUES
    xmat = np.dot(reclat,kpoints)
    lengthx = len(xmat[0,:])
    xvalues = np.zeros(lengthx)
    total = 0
    for i in range(1,lengthx):
	distance = LA.norm(xmat[:,i]-xmat[:,i-1])
	total+=distance
	xvalues[i] = total

    #Inputs (Minimums and Maximums)
    xmax = max(xvalues)
    ymin = fermilevel-10
    ymax = fermilevel+10
    #ymin = np.amin(ymat)
    #ymax = np.amax(ymat)

    #Define Full Figure
    width=14
    height=7.5
    fig1 = plt.figure(figsize=(width,height),dpi=200)
    gs = gridspec.GridSpec(2,7)
    matplotlib.rc('xtick',labelsize=20)
    matplotlib.rc('ytick',labelsize=20)

    #Define Band Plot Location and Attributions
    ax1 = plt.subplot(gs[:,0:5])
    ax1.tick_params(axis='x',which='both',bottom='on', top='off',labelbottom='off')
    ax1.axis([0,xmax,ymin,ymax])
    ax1.set_ylabel('Energy (eV)',size=20)

    #SUMMARIZE CHARACTERS
    characters_new = summarize_characters(characters)

    #Plot Electronic Bands
    if ymat.shape[0] == 2: #SPIN POLARIZED
	for n in range(len(ymat[0,0,:])):
	    nup_plot = rgbline(ax1, xvalues, ymat[0,:,n], characters_new[:,n,1], characters_new[:,n,2], characters_new[:,n,3], alpha=1.)
	    ndown_plot = rgbline(ax1, xvalues, ymat[1,:,n], characters_new[:,n,1], characters_new[:,n,2], characters_new[:,n,3], alpha=1.)		    
    else:
	for n in range(len(ymat[0,0,:])):	    
	    nplot = rgbline(ax1,xvalues,ymat[0,:,n],characters_new[:,n,1],characters_new[:,n,2],characters_new[:,n,3],alpha=1.)	

    #Set Dimensions of Image Canvas
    start=0.06
    finish=0.65
    scale = finish-start
    plt.subplots_adjust(left=start+0.01,right=0.9,bottom=0.1,top=0.92)
    
    #Generate First and Last Label
    plt.figtext(start,0.05,labels[0],size=20)
    plt.figtext(start+scale,0.05,labels[len(labels)-1],size=20)
    
    #Set Remaining Partitions and Labels
    x_spec = np.dot(reclat,kpt_spec)
    distance=0
    for m in range(len(x_spec[0,:])-1):
	interval = LA.norm(x_spec[:,m+1]-x_spec[:,m])
	distance += interval
	ax1.vlines(distance,ymin,ymax,colors='black',linewidth=2)
	plt.figtext(start+scale*distance/max(xvalues),0.05,labels[m+1],size='20')

    #Plot Fermi Level
    fermi_plot = ax1.hlines(fermilevel,0,xmax,colors='orange',linestyles=u'dashed',linewidth=2,label='Fermi Level')
    
    #LEGEND
    splot = plt.hlines(-2000,0,xmax,color='k',label='s-band')
    pplot = plt.hlines(-2000,0,xmax,color='r',label='p-band')
    dplot = plt.hlines(-2000,0,xmax,color='g',label='d-band')
    fplot = plt.hlines(-2000,0,xmax,color='b',label='f-band')    

    #Define DOS Plot Location and Attributes
    ax2 = plt.subplot(gs[:,5:])
    ax2.tick_params(axis='x',which='both',bottom='on', top='on',labelbottom='off')
    ax2.axes.get_yaxis().set_ticklabels([])

    #Plot DOS
    summarized_dosy = summarize_dos_characters(dosy,total_dosy)
    
    print total_dosy.shape[0]
    if total_dosy.shape[0]==2:
	total_dosy[0,:] = total_dosy[0,:]+total_dosy[1,:]
    xmax2 = np.amax(total_dosy)
    ax2.axis([0,xmax2,ymin,ymax])    

    ax2.plot(summarized_dosy[:,0],dosx[0,:],color='k',linewidth=2,label="s-dos")
    ax2.plot(summarized_dosy[:,1],dosx[0,:],color='r',linewidth=2,label="p-dos")
    ax2.plot(summarized_dosy[:,2],dosx[0,:],color='g',linewidth=2,label="d-dos")
    ax2.plot(summarized_dosy[:,3],dosx[0,:],color='b',linewidth=2,label="f-dos")
    ax2.fill_betweenx(dosx[0,:],0,total_dosy[0,:],color=(0.7,0.7,0.7),facecolor=(0.7,0.7,0.7))
    ax2.plot(total_dosy[0,:],dosx[0,:],color=(0.6,0.6,0.6))

    #Plot Fermi Level 
    ax2.hlines(fermilevel,0,xmax2,colors='orange',linestyles=u'dashed',linewidth=2)

    #Plot Legend
    #Legend
    plt.legend(handles=[splot,pplot,dplot,fplot],bbox_to_anchor=(1.4,1.01))     

    #Plot Save
    plt.savefig(dir_save+"/Orbital_ElectronicBSandDOS.pdf")
    return 

#===========================  PHONON BS AND DOS PLOT  ====================================   
def phonbanddosplot(qpoints,bands,dosx,dosy,reclat,kspec,labels,dir_save):
	
	#Inputs - Convert from eV to meV
	hbar = 0.004135
	ymin = 0
	ymax = 1000*hbar*max(dosx)
	xmax2 = max(dosy)*1.2
	
	#GEN X-VALUES
	xmat = np.dot(reclat,qpoints)
	lengthx = len(xmat[0,:])
	xvals = np.zeros(lengthx)
	total = 0
	for i in range(1,lengthx):
	    distance = LA.norm(xmat[:,i]-xmat[:,i-1])
	    total+=distance
	    xvals[i] = total	
	    
	xmax1 = max(xvals)

	#Generate Phonon Dispersion Plot
	width=14
	height=7.5
	fig = plt.figure(figsize=(width,height),dpi=200)
	gs = gridspec.GridSpec(2,7)
	matplotlib.rc('xtick',labelsize=20)
	matplotlib.rc('ytick',labelsize=20)
	ax1 = plt.subplot(gs[:,0:5])
	ax1.tick_params(axis='x',which='both',bottom='on', top='off',labelbottom='off')
	for i in range(bands.shape[0]):
		if (i<3):
			band_plot_ac,=ax1.plot(xvals,1000*bands[i,:],color='k',linewidth=3,label='Acoustic')
		else:
			band_plot_op,=ax1.plot(xvals,1000*bands[i,:],color='gray',linestyle='--',linewidth=3,label='Optical')
    
	ax1.axis([0,xmax1,ymin,ymax])
	ax1.set_ylabel('Energy (meV)',size=20)
	plt.title('Phonon Dispersion and DOS',size=35,color='black')

	#Set Dimensions of Image Canvas
	start=0.065
	finish=0.65
	scale = (finish-start)
	plt.subplots_adjust(left=start+0.01,right=0.9,bottom=0.1,top=0.92)
	
	#Generate First and Last Label
	plt.figtext(start,0.05,labels[0],size=20)
	plt.figtext(start+scale,0.05,labels[len(labels)-1],size=20)
	
	#Set Remaining Partitions and Labels
	x_spec = np.dot(reclat,kspec)
	distance=0
	for m in range(len(x_spec[0,:])-1):
	    interval = LA.norm(x_spec[:,m+1]-x_spec[:,m])
	    distance += interval
	    ax1.vlines(distance,ymin,ymax,colors='black',linewidth=2)
	    plt.figtext(start+scale*distance/max(xvals),0.05,labels[m+1],size='20')	

	#Generate DOS Plot
	ax2 = plt.subplot(gs[:,5:])
	ax2.plot(dosy,1000*hbar*dosx,color='darkred',linewidth=3)
	ax2.tick_params(axis='x',which='both',bottom='on', top='on',labelbottom='off')
	ax2.axes.get_yaxis().set_ticklabels([])
	ax2.axis([0,xmax2,ymin,ymax])

	#Legend
	if (i<3):
		plt.legend(handles=[band_plot_ac],bbox_to_anchor=(1.4,1.01))
	else:
		plt.legend(handles=[band_plot_ac,band_plot_op],bbox_to_anchor=(1.4,1.01))
        
	plt.savefig(dir_save+"/PhononBSandDOS.png")
	return

#=========================== PHONON BS PLOT  =============================================
def phonbandplot(qpoints,bands,reclat,kspec,labels,dir_save):
	
	#INPUTS - Convert from eV to meV
	ymin = 1000*np.amin(bands)
	ymax = 1000*np.amax(bands)
	
	#GEN X-VALUES
	xmat = np.dot(reclat,qpoints)
	lengthx = len(xmat[0,:])
	xvals = np.zeros(lengthx)
	total = 0
	for i in range(1,lengthx):
	    distance = LA.norm(xmat[:,i]-xmat[:,i-1])
	    total+=distance
	    xvals[i] = total	
    
	xmax = max(xvals)

	#Generate Phonon Dispersion Plot
	width=10
	height=7
	fig = plt.figure(figsize=(width,height),dpi=200)
	matplotlib.rc('xtick',labelsize=20)
	matplotlib.rc('ytick',labelsize=20)
	ax = plt.subplot()
	ax.tick_params(axis='x',which='both',bottom='on', top='off',labelbottom='off')
	for i in range(len(bands)):
		if (i<3):
			band_plot_ac,=ax.plot(xvals,1000*bands[i,:],color='k',linewidth=3,label='Acoustic')
		else:
			band_plot_op,=ax.plot(xvals,1000*bands[i,:],color='gray',linestyle='--',linewidth=3,label='Optical')
    
	ax.axis([0,xmax,ymin,ymax])
	ax.set_ylabel('Energy (meV)',size=20)
	plt.title('Phonon Dispersion',size=35,color='black')

	#Set Dimensions of Image Canvas
	start=0.08
	finish=0.96
	scale=(finish-start)
	plt.subplots_adjust(left=start+0.01,right=finish,bottom=0.1,top=0.92)
	
	#Generate First and Last Label
	plt.figtext(start,0.05,labels[0],size=20)
	plt.figtext(start+scale,0.05,labels[len(labels)-1],size=20)
    
	#Set Remaining Partitions and Labels
	x_spec = np.dot(reclat,kspec)
	distance=0
	for m in range(len(x_spec[0,:])-1):
	    interval = LA.norm(x_spec[:,m+1]-x_spec[:,m])
	    distance += interval
	    ax.vlines(distance,ymin,ymax,colors='black',linewidth=2)
	    plt.figtext(start+scale*distance/max(xvals),0.05,labels[m+1],size='20')	
		
	#Legend
	if (i<3):
		plt.legend(handles=[band_plot_ac])
	else:
		plt.legend(handles=[band_plot_ac,band_plot_op])
		
	plt.savefig(dir_save+"/PhononBS.png")
	return

#=============================  PHONON DOS PLOT  ==========================================
def phondosplot(dosx,dosy,dir_save):
    
    #Define Full Figure
    fig = plt.figure(figsize=(10,7),dpi=200)
    matplotlib.rc('xtick',labelsize=20)
    matplotlib.rc('ytick',labelsize=20)

    #Min and Max
    hbar = 0.004135
    xmin=1000*hbar*min(dosx)
    xmax=1000*hbar*max(dosx)
    ymin=np.amin(dosy)
    ymax=np.amax(dosy)

    #Define DOS Plot Location and Attributes
    ax1 = plt.subplot()
    ax1.tick_params(axis='x',which='both',bottom='on', top='on',labelbottom='on')
    ax1.axes.get_yaxis().set_ticklabels([])
    ax1.axis([xmin,xmax,ymin,ymax])
    ax1.grid(True)
    ax1.set_ylabel('Amplitude',size=20)
    ax1.set_xlabel('Energy (meV)',size=20)
    plt.title('Phonon DOS',size=35,color='black')

    #Plot DOS
    ax1.plot(1000*hbar*dosx,dosy,color='darkred',linewidth=3,label="Phonon DOS")
    plt.legend()
       
    plt.savefig(dir_save+'/PhononDOS.png')
    return

def conv(xvalues,yvalues,time,ymin,ymax,tmin,tmax,parameter,label,dir_save):
    
    #Define Figure Parameters
    fig = plt.figure(figsize=(9,6.5),dpi=200)
    matplotlib.rc('xtick',labelsize=20)
    matplotlib.rc('ytick',labelsize=20)
    
    #Generate Time Plot
    fillarray=np.zeros(len(time))
    ax1=plt.subplot()
    time2, = ax1.plot(xvalues, time, '-o',color='orange',linewidth=2,label="TIME")
    ax1.fill_between(xvalues, fillarray,time,color='lightgray')
    ax1.set_ylim(tmin, tmax); #Min/Max
    ax1.set_xlabel(parameter,size=15)
    ax1.set_ylabel('Total CPU Time (secs)',size=15)
    ax1.set_title('Convergence Plot',size=20)
    
    #Generate Converge Plot
    ax2 = ax1.twinx()
    value, = ax2.plot(xvalues, yvalues,'-o',color='k',linewidth=2,label="VALUE")
    ax2.grid(True)
    ax2.set_ylim(ymin,ymax)
    ax2.set_ylabel(label,size=15)
    plt.legend(handles=[value,time2])
    
    plt.savefig(dir_save+'/'+parameter+'.png')

#===============================       ===================================================
#===============================  END  ===================================================
#===============================       ===================================================