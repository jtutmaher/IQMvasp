## JAKE A TUTMAHER
## JOHNS HOPKINS UNIVERSITY
## MCQUEEN LABORATORY
## THE INSTITUTE FOR QUANTUM MATTER
## DEPARTMENT OF PHYSICS, DEPARTMENT OF CHEMISTRY, DEPARTMENT OF MATERIALS SCIENCE AND ENGINEERING
##
## CONTACT: jtutmah1@jhu.edu
###################################################################################################

""" Convert VASP EIGENVAL file to bxsf file for XCRYSDEN"""

#=================  MODULES  ========================
import os
import numpy as np
import glob
import IQM.VASPread as VASPread
import utils.inputs as inputs

#=================  Directories  ====================
#RETURN CURRENT / PARENT DIRECTORY
dir1 = inputs.get_current_directory()

# !!!!!!!!!! SET THIS VALUE !!!!!!!!!!!!!!!!!!!!!!!
dir_ele = dir1+'STATIC/'

#INCAR FILE FOR SOC AND SPINPOL
incarfile = VASPread.incar(dir_ele)
spinpol=incarfile.spin()
soc=incarfile.soc()

#INITIATE EIGENVAL AND OUTCAR FILES
eigenvalfile = VASPread.eigenval(dir_ele,SpinPol=spinpol,SOC=soc)
outfile = VASPread.outcar(dir_ele)

#PULL RELEVANT DATA USING VASPread
kpoints = eigenvalfile.kpoints()
nband = eigenvalfile.nband()

reclat = outfile.reclatvec()
name = outfile.compound()
nkpts = outfile.nkpts()
fermilevel = outfile.fermilevel()

#DETERMINE KPT GRID FROM NKPTS
singlek = int(round(nkpts**(1./3)))
ksamp = [singlek,singlek,singlek]

#GET BAND ENERGIES
energies = eigenvalfile.energy()

#REFORMAT ENERGIES AND GENERATE KMESH
new_energy = inputs.reformat_energy(energies)
kmat,energymat = inputs.meshgrid(kpoints,new_energy)

#GENERATE BXSF FILE
if os.path.exists(dir1+name+'.bxsf'):
    os.remove(dir1+name+'.bxsf')

#GENERATE BXSF HEADER    
bxsf = open(name+'.bxsf', mode='a')
bxsf.write('BEGIN_INFO \n  #Created by Jake Tutmaher and Guy Marcus\n  #Not for outside distribution\n')
bxsf.write('  Fermi_Energy: '+str(fermilevel)+'\n')
bxsf.write('END_INFO\n\nBEGIN_BLOCK_BANDGRID_3D\n  Pointless_Line\n  BEGIN_BANDGRID_3D\n')
bxsf.write('    '+str(nband)+'\n')
bxsf.write('    '+str(ksamp[0])+' '+str(ksamp[1])+' '+str(ksamp[2])+'\n    0.0 0.0 0.0\n')
bxsf.write('    '+str(reclat[0,0])+' '+str(reclat[0,1])+' '+str(reclat[0,2])+'\n')
bxsf.write('    '+str(reclat[1,0])+' '+str(reclat[1,1])+' '+str(reclat[1,2])+'\n')
bxsf.write('    '+str(reclat[2,0])+' '+str(reclat[2,1])+' '+str(reclat[2,2])+'\n')

#GENERATE ENERGIES IN ROW MAJOR
for n in range(nband):
    bxsf.write('\n  BAND: '+str(n+1)+'\n')
    for m in range(ksamp[0]):
        for p in range(ksamp[1]):
            for q in range(ksamp[2]):
                if abs(energymat[q,p,m,n])<10:
                    energy2 = format(energymat[m,p,q,n],'.5f')
                else:
                    energy2 = format(energymat[m,p,q,n],'.4f')
                bxsf.write('    '+str(energy2))
            bxsf.write('\n')
        if m==(ksamp[0]-1):
            continue
        else:
            bxsf.write('\n\n')
bxsf.write('  END_BANDGRID_3D\n')
bxsf.write('END_BLOCK_BANDGRID_3D')
bxsf.close()
