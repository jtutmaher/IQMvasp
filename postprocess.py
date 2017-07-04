#Jacob A. Tutmaher
#July 15, 2015
#Johns Hopkins University, Department of Physics
#
#Main Script which calls several functions from Modules IQMvasp and IQMplots to generate 
#BS and DOS plots. IQMvasp reads and outputs Electronic and Phonon data from VASP and 
#Phonopy files into numpy arrays. IQMplots talks numpy arrays, floats, and directories to 
#generate and save plots. There are no optional arguments for these functions, so follow
#the help directions exactly. A summary file with potentially relevant run information is 
#also generated at the end of this document.
#
#Contact: jtutmah1@jhu.edu

#==================================  MODULES  ============================================

from IQM import VASPread
from IQM import PHONOPYread
from IQM import plots
import numpy as np
import os
#import IQM.vasp as vasp

#=================================  DIRECTORY  ===========================================

#Returns Current Directory
dir1 = os.popen("pwd").readlines()[0]
length=len(dir1)
dir1 = dir1[:length-1]

counter=0;

#Check that filepaths exist
if not os.path.exists(dir1+'/PATH/'):
	print('Could Not Find PATH Directory')
	counter=counter+1
if not os.path.exists(dir1+'/STATIC/'):
	print('Could Not Find STATIC Directory')
	counter=counter+1

#Define filepaths	
dir_ele = dir1+'/PATH/'
dir_eledos = dir1+'/STATIC/'
dir_phon = dir1+'/DFPT/'

#Save directories
if not os.path.exists(dir1+'/PostProcess/') and not counter>1:
	os.mkdir(dir1+'/PostProcess/')	

#WRITE SUMMARY FILE
dir_save = dir1+'/PostProcess/'
summary=open(dir_save+"/Summary.txt",'w')

#=================================  ELECTRONIC ===========================================

if os.path.exists(dir_ele):
	
	#INITIAL CONDITIONS FOR PLOTTING ENVIRONMENT
	if os.path.exists(dir_ele+'PROCAR') and os.path.exists(dir_eledos+'PROCAR'):
		procar=True
	else:
		procar=False
	
	incarfile = VASPread.incar(dir_eledos)
	spinpol=incarfile.spin()
	soc=incarfile.soc()
	
	#INITIALIZE FILES
	outfile = VASPread.outcar(dir_ele)
	eigenfile = VASPread.eigenval(dir_ele,SpinPol=spinpol,SOC=soc)
	kfile = VASPread.kpoints(dir_ele)
	dosfile = VASPread.doscar(dir_eledos,SpinPol=spinpol,SOC=soc)
	
	#DIRECT AND RECIPROCAL LATTICE INFORMATION	
	kpoints = eigenfile.kpoints()
	dirlat = outfile.dirlatvec()
	reclat = outfile.reclatvec()
	
	#FERMI LEVEL
	fermilevel=outfile.fermilevel()
	
	#ENERGIES
	energies = eigenfile.energy()
	
	#BAND PATH LABELS
	kspec = kfile.kpath()
	labels = kfile.labels()
	
	#DOS
	dosx = dosfile.energy()
	dosy = dosfile.dos()

	#ELECTRONIC PLOTS
	if len(kpoints[0,:])==len(energies[0,:,0]):
		#STANDARD BAND PLOTS
		elestr='Dimensions Agree, Generated Electronic Plots'
		plots.elebanddosplot(energies,reclat,kpoints,fermilevel,kspec,labels,dosx,dosy,dir_save)
		plots.elebandplot(energies,reclat,kpoints,fermilevel,kspec,labels,dir_save)
		plots.eledosplot(dosx,dosy,fermilevel,dir_save)
		
		if procar==True:
			
			print "PARSING DOSCAR AND PROCAR FILE"
			
			#ORBITAL INFORMATION
			profile = VASPread.procar(dir_ele,SpinPol=spinpol,SOC=soc)
			oenergy,odos = dosfile.odos()
			kpoints_orbital = profile.kpoints()
			energies_orbital = profile.energies()
			characters = profile.character()
			
			#ORBITAL PLOTS
			plots.orbital_dosplot(oenergy,odos,dosy,fermilevel,dir_save)
			plots.orbital_elebandplot(energies_orbital,reclat,kpoints_orbital,fermilevel,kspec,labels,characters,dir_save)
			plots.orbital_elebanddosplot(energies_orbital,reclat,kpoints_orbital,fermilevel,kspec,labels,characters,oenergy,odos,dosy,dir_save)
	else:
		elestr='Dimenstions Do Not Agree, Did Not Generate Electronic Plots'
	
	summary.write('----ELECTRONIC ANALYSIS----'+'\n')	
	summary.write('The Direct Lattice Is'+'\n'+str(dirlat)+'\n')
	summary.write('The Reciporcal Lattice Is'+'\n'+str(reclat)+'\n')
	summary.write('The Brillouin Zone Path Is'+'\n'+str(labels)+'\n'+str(kspec)+'\n')
	summary.write(elestr+'\n \n')	
else:
	pass
#===================================  PHONONS  ===========================================

if os.path.exists(dir_phon):
	##DISPERSION
	#phonx, phonband, interval2 = vasp.phonbands(dir_phon)
	
	#TEST
	phonfile = PHONOPYread.bands(dir_phon)
	qbands = phonfile.bands()
	qpoints = phonfile.qpoints()
	#labels=['T','T','T','T','T','T','T']

	#DOS
	meshfile = PHONOPYread.mesh(dir_phon)
	qdosx,qdosy = meshfile.dos()
	#phondosx, phondosy = vasp.phondos(dir_phon)
	
	#SPECIAL POINTS
	qspec = kspec

	#PHONON PLOTS
	if qbands.shape[1]==qpoints.shape[1] and qdosx.shape==qdosy.shape:
		plots.phonbanddosplot(qpoints,qbands,qdosx,qdosy,reclat,qspec,labels,dir_save)
		plots.phonbandplot(qpoints,qbands,reclat,qspec,labels,dir_save)
		plots.phondosplot(qdosx,qdosy,dir_save)
		phonstr='Dimensions Agree, Generated Phonon Plots'
	else:
		phonstr='Dimensions Do Not Agree, Did Not Generate Phonon Plots'
	
	summary.write('----PHONON ANALYSIS----'+'\n')
	summary.write(phonstr+'\n')
	summary.close()	
else:
	pass


#====================================       ==============================================
#====================================  END  ==============================================
#====================================       ==============================================