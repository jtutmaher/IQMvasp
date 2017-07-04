## JAKE A TUTMAHER
## JOHNS HOPKINS UNIVERSITY
## MCQUEEN LABORATORY
## THE INSTITUTE FOR QUANTUM MATTER
## DEPARTMENT OF PHYSICS, DEPARTMENT OF CHEMISTRY, DEPARTMENT OF MATERIALS SCIENCE AND ENGINEERING
##
## CONTACT: jtutmah1@jhu.edu
###################################################################################################

"""PHONOPY specific readers - generating the specified information from the raw yaml files. For 
   further information, please see www.moreinfo.com. NOTE: Requires a python yaml reader and 
   numpy.
   
   ---
   
   qpoints: class which parses information in Phonopy qpoints.yaml file and returns it as a numpy
            array - or float depending. Requires the use of a band.conf (with QPOINTS=.TRUE.) 
            input file beforehand - contains phonon energy/frequency information. Returns energies
            (in eV) - NOT frequency (in THz).
   mesh:    class which parses information in Phonopy mesh.yaml file and returns it as a numpy array.
            Requires use of a mesh.conf input file beforehand - contains DOS information.
   band:    class which parses information in Phonopy band.yaml file and returns it as a numpy array
            - or float depending. Requires the use of a band.conf input file with q-path defined
            beforehand. Returns energies in eV - NOT frequency (in THz).
"""
####################################################################################################

#===========================    IMPORT SPECIFIC PACKAGES   =================================
import numpy as np
import yaml
import os

#============================    DEFINE QPOINTS CLASS   =====================================
class qpoints:
    """ Qpoints class read the qpoints.yaml file contained in a DFPT directory."""
    def __init__(self,phon_dir):
        self.phon_dir = phon_dir
        with open(self.phon_dir+'qpoints.yaml','r') as f:
            qpoints = yaml.load(f)
        self.qpoints = qpoints
        return
    
    def reclat(self):
        rfile = np.zeros((3,len(self.qpoints["reciprocal_lattice"])))
        for i in range(len(self.qpoints["reciprocal_lattice"])):
            rfile[:,i] = self.qpoints["reciprocal_lattice"][i]
        return rfile
    
    def kpoints(self):       
        qfile = np.zeros((3,len(self.qpoints["phonon"])))
        for i in range(len(self.qpoints["phonon"])):
            qfile[:,i] = self.qpoints["phonon"][i]["q-position"]
        return qfile
    
    def bands(self):
        hbar = 0.004135 #eV/THz       
        bandfile = np.zeros((len(self.qpoints["phonon"]),len(self.qpoints["phonon"][0]["band"])))
        for i in range(len(self.qpoints["phonon"])):
            for j in range(len(self.qpoints["phonon"][0]["band"])):
                bandfile[i,j] = self.qpoints["phonon"][i]["band"][j]["frequency"]
        bandfile = hbar*bandfile
        return bandfile
    
    def natoms(self):
        natoms = self.qpoints["natom"]
        return natoms
    
    def nbands(self):
        nbands = len(self.qpoints["phonon"][0]["band"])
        return nbands

#=============================    DEFINE MESH CLASS   ======================================
    
class mesh:
    """ The mesh class reads the total_dos.dat and dos files. Could be improved
        to actually read mesh.yaml."""
    def __init__(self,dos_dir):
        self.dos_dir = dos_dir
        if not os.path.exists(self.dos_dir+'total_dos.dat'):
            os.chdir(self.dos_dir)
            os.system('phonopy -c POSCAR-unitcell --dos mesh.conf')
        else:
            pass
        dosstr = open(self.dos_dir+'total_dos.dat','r')
        self.dosstr = dosstr.read()
        dosstr.close()
        return
        
    def dos(self):
        dosstr2 = self.dosstr
        dosstr2 = dosstr2.split('\n')[1:]
        energy = []
        dos = []
        for i in range(len(dosstr2)-1):
            energy = np.append(energy,float(dosstr2[i].split()[0]))
            dos = np.append(dos,float(dosstr2[i].split()[1]))
        hbar = 0.004135 #eV/THz
        #energy = hbar*energy
        return energy,dos

#============================    DEFINE BANDS CLASS   ======================================
    
class bands:
    """ Class to read band.yaml file in the DFPT directory."""
    def __init__(self,band_dir):
        self.band_dir = band_dir
        with open(self.band_dir+'band.yaml','r') as f:
            band = yaml.load(f)
        self.band = band
        return
    
    def reclat(self):
        rfile = np.zeros((3,len(self.band["reciprocal_lattice"])))
        for i in range(len(self.band["reciprocal_lattice"])):
            rfile[:,i] = self.band["reciprocal_lattice"][i]
        return rfile
    
    def qpoints(self):       
        qfile = np.zeros((3,len(self.band["phonon"])))
        for i in range(len(self.band["phonon"])):
            qfile[:,i] = self.band["phonon"][i]["q-position"]
        return qfile
    
    def distance(self):
        dfile = np.zeros((1,len(self.band["phonon"])))
        for i in range(dfile.shape[1]):
            dfile[0,i] = self.band["phonon"][i]["distance"]
        return dfile
    
    def bands(self):
        hbar = 0.004135 #eV/THz 
        bandfile = np.zeros((len(self.band["phonon"][0]["band"]),len(self.band["phonon"])))
        for i in range(len(self.band["phonon"][0]["band"])):
            for j in range(len(self.band["phonon"])):
                bandfile[i,j] = self.band["phonon"][j]["band"][i]["frequency"]
        bandfile = hbar*bandfile
        return bandfile
    
    def natoms(self):
        natoms = self.band["natom"]
        return natoms
    
    def nbands(self):
        nbands = len(self.band["phonon"][0]["band"])
        return nbands    
    
    
        
        
        