## Jake A Tutmaher
## JOHNS HOPKINS UNIVERSITY
## MCQUEEN LABORATORY
## DEPARTMENT OF PHYSICS, CHEMISTRY, MATERIALS SCIENCE
##
##
## July 16, 2016 
##
######################################################################
"""Generate Gamma Centered Meshes for Fermi Plots in Mayavi"""

#IMPORT
import os
import numpy as np
import utils.inputs as inputs
import IQM.VASPread as VASPread

#GET CURRENT DIRECTORY
dir1 = inputs.get_current_directory()

# !!!!!!!!!! SET THIS VALUE !!!!!!!!!!!!!!!!!!!!!!!
dir_ele = dir1+'STATIC/'

#INCAR FILE FOR SOC AND SPINPOL
incarfile = VASPread.incar(dir_ele)
spinpol=incarfile.spin()
soc=incarfile.soc()

#INITIATE RELEVANT OUTPUT FILES
outfile = VASPread.outcar(dir_ele)
eigenfile = VASPread.eigenval(dir_ele,SpinPol=spinpol,SOC=soc)

#GET DATA USING VASPREAD
energy = eigenfile.energy()
kpts = eigenfile.kpoints()
fermilevel = outfile.fermilevel()

#REFORMAT DATA
new_energy = inputs.reformat_energy(energy)
kmat,final_energy = inputs.meshgrid(kpts,new_energy)

#SAVE
np.savez(dir1+"surface.npz",kpoints=kmat,energy=final_energy,fermilevel=fermilevel)