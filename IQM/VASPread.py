## JAKE A TUTMAHER
## JOHNS HOPKINS UNIVERSITY
## MCQUEEN LABORATORY   
## THE INSTITUTE FOR QUANTUM MATTER
## DEPARTMENT OF PHYSICS, DEPARTMENT OF CHEMISTRY, DEPARTMENT OF MATERIALS SCIENCE AND ENGINEERING
##
## CONTACT: jtutmah1@jhu.edu
##################################################################################################

"""VASP specific readers - generating the specified information from the raw output files. 
   For further information, please see www.moreinfo.com.
   
   ---
   
   eigenval:   class which parses VASP eigenval file to return kpoints and energies (in eV)
               as numpy arrays or floats - depending. Flag for Spin Polarized vs Non- Spin 
               Polarized (and SOC) depending - as VASP format changes.
   outcar:     class which parses VASP outcar file and returns system information, such as
               fermilevel, direct lattice vectors, reciprocal lattice vectors and compound.
               Returns floats, strings, or numpy arrays depending.
   doscar:     class which parses VASP doscar file and returns associated energies and DOS.
               Contains methods that can also parse orbital character information in the 
               event that LORBIT is set in VASP INCAR. Capable of handle both Spin Polarized
               and Non-Sping Polarized (plus SOC) cases as the output format changes between
               cases.
   procar:     class which parses VASP procar file. This file only exists if LORBIT is set in
               VASP INCAR file. For SOC case, it returns character information for all four
               spinors. For the Spin Polarized case it returns character information for both
               spin-up and spin-down cases. Else - it returns simply the total character at
               each k-point and energy level. Returns information in the form of numpy array.
   qscript:    Returns information contained in a bash queing script. Design for SLURM queue
               system. May not be applicable for all users.
"""

###############################################################################################

#===========================    IMPORT SPECIFIC PACKAGES   =================================
import numpy as np
import os
import re

#===========================    DEFINE EIGENVAL CLASS   ====================================
class eigenval:
    """ Eigenval class reads stated information from the EIGENVAL file in given directory."""
    
    def __init__(self,dir_eig,SpinPol=True,SOC=True):
        #DEFINE EIGENVAL FILE AND OPEN AS STRING
        self.dir1 = dir_eig
        eigenstr = open(self.dir1+'EIGENVAL','r')
        eigenstr2 = eigenstr.read().split('\n \n')
        eigenstr.close()
        
        #DEFINE EIGENVAL HEADER AND BODY
        self.header = eigenstr2[0]
        self.eigenval = eigenstr2[1:]
        fsteps = self.header.split('\n')[5]
        self.nbands = int(fsteps.split()[2])
        
        #DEFINE IF SPINPOL OR SOC LAYOUT
        self.spinpol = SpinPol
        self.soc = SOC
        return
        
    def kpoints(self):
        #READ KPOINTS FROM EIGENVAL (3D VECTOR)
        kpts = np.zeros([3,len(self.eigenval)])
        for i in range(len(self.eigenval)):
            kpts_prelim = self.eigenval[i].split('\n')
            kpts_str = kpts_prelim[0]
            kx = kpts_str.split()[0]
            ky = kpts_str.split()[1]
            kz = kpts_str.split()[2]
            kvec = [float(kx),float(ky),float(kz)]
            kpts[:,i]=kvec   
        return kpts
    
    def nband(self):
        return self.nbands
    
    def nvalence(self):
        #READ VALANCE ELECTRONS FROM EIGENVAL HEADER
        valencestr = self.header.split('\n')[5]
        valence = valencestr.split()[0]
        nvalence = int(valence)
        return nvalence
    
    def energy(self): 
        if self.spinpol==True and self.soc==False:
            #READ ENERGY FOR SPIN POLARIZED (SPIN UP FIRST, THEN SPIN DOWN)
            energy_array = np.zeros((2,len(self.eigenval),self.nbands))
            for i in range(len(self.eigenval)):
                line = self.eigenval[i].split('\n')
                e_str = line[1:]
                for j in range(self.nbands):
                    energies = e_str[j].split()
                    energy_array[0,i,j] = float(energies[1])
                    energy_array[1,i,j] = float(energies[2])
        else:
            #READ ENERGY FOR NON-SPIN POLARIZED - IN ORDER
            energy_array = np.zeros((1,len(self.eigenval),self.nbands))
            for i in range(len(self.eigenval)):
                line = self.eigenval[i].split('\n')
                e_str = line[1:]
                for j in range(self.nbands):
                    energies = e_str[j].split()
                    energy_array[0,i,j] = float(energies[1])                
        
        return energy_array

#===========================    DEFINE OUTCAR CLASS   ====================================       
class outcar:
    """ Class to read parameters from OUTCAR file."""
    
    def __init__(self,dir_outcar):
        #READ OUTCAR FILE AS STRING
        self.dir_outcar = dir_outcar
        outstr = open(self.dir_outcar+'OUTCAR','r')
        outstr2 = outstr.read()
        outstr.close()
        self.outstr = outstr2
        return
        
    def latticeconst(self):
        #FIND AND READ LATTICE VECTORS FROM OUTCAR AS STRING
        start = self.outstr.find("ALAT")
        stop = self.outstr.find("Lattice vectors:")
        latticeconst = float(self.outstr[start:stop].split()[2])  
        return latticeconst
    
    def system(self):
        #READ SYSTEM (COMPOUND NAME) AS STRING FROM OUTCAR
        start = self.outstr.find("SYSTEM =")
        stop = self.outstr.find("POSCAR =")
        name = self.outstr[start:stop].split()[2]
        return name
    
    def dirlatvec(self):
        #READ DIRECT LATTICE VECTORS FROM OUTCAR - RETURN ARRAY OF FLOATS
        start = self.outstr.find("direct lattice vectors")
        stop = self.outstr.find("length of vectors")
        latticestr = self.outstr[start:stop].split("\n")
        dirlat=np.zeros((3,3))
        for i in range(1,4):
            latticestr2 = latticestr[i].split()
            dirlat[0,i-1] = float(latticestr2[0])
            dirlat[1,i-1] = float(latticestr2[1])
            dirlat[2,i-1] = float(latticestr2[2])
        return dirlat
    
    def reclatvec(self):
        #READ RECIPROCAL LATTICE VECTORS FROM OUTCAR - RETURN ARRAY OF FLOATS
        start = self.outstr.find("direct lattice vectors")
        stop = self.outstr.find("length of vectors")
        latticestr = self.outstr[start:stop].split("\n")
        reclat=np.zeros((3,3))
        for i in range(1,4):
            latticestr2 = latticestr[i].split()
            reclat[0,i-1] = float(latticestr2[3])
            reclat[1,i-1] = float(latticestr2[4])
            reclat[2,i-1] = float(latticestr2[5]) 
        return reclat
    
    def fermilevel(self):
        #READ FERMILEVEL FROM OUTCAR - SOMETIMES OUTCAR DOESNT CONTAIN THIS AND WILL THROW AN ERROR
        #CHECK OUTCAR FILE FIRST BEFORE ADJUSTING THIS METHOD
        start = self.outstr.find("E-fermi")
        stop = self.outstr.find("XC(G=0)")
        fermifile = self.outstr[start:stop].split()
        endpoint=len(fermifile)
        fermilevel=float(fermifile[endpoint-1])
        return fermilevel
    
    def compound(self):
        #SAME AS SYSTEM
        start = self.outstr.find("SYSTEM =")
        stop = self.outstr.find("POSCAR =")
        compound = self.outstr[start:stop].split()[2]
        return compound
    
    def nions(self):
        #GET NUMBER OF ATOMS IN UNIT CELL
        start = self.outstr.find("NIONS =")
        stop = self.outstr.find("non local maximal ")
        nions = self.outstr[start:stop].split()[2]
        return float(nions)
    
    def nkpts(self):
        #GET TOTAL NUMBER OF KPOINTS
        start = self.outstr.find("NKPTS =")
        stop = self.outstr.find("k-points in BZ")
        nkpts = float(self.outstr[start:stop].split()[2])
        return nkpts                
    
    def time(self):
        #GET TOTAL TIME NEEDED TO COMPUTE VASP JOB
        try:
            location = self.outstr.find('Total CPU time used (sec):')
            finalstring = self.outstr[location:].split('\n')[0]
            finalstring = finalstring.split()[5]
            time = float(finalstring)/3600 #Convert to Hours
            return time
        except IndexError:
            time = 0
            return time
    
#===========================    DEFINE DOSCAR CLASS   ====================================
class doscar:
    """ Class to read doscar file and pull DOS/Energy values. It may be easier to read 
        this info correctly from vasprun.xml."""
    
    def __init__(self,dir_doscar,SpinPol=True,SOC=True):
        #READ IN DOSCAR AS STRING
        self.dir_doscar = dir_doscar
        dosfile = open(self.dir_doscar+'DOSCAR','r')
        dosstr = dosfile.read()
        dosfile.close()
        self.dosstr=dosstr
        #GET LENGTH OF DOSSTR
        line = self.dosstr.split('\n')[5]
        self.length = int(line.split()[2])
        self.soc=SOC
        self.spinpol=SpinPol
        return
        
    def dos(self):
        #READ DOS AS ARRAY OF FLOATS
        dosstr = self.dosstr.split('\n')
        count=0
        if self.soc==False and self.spinpol==True:
            #FIRST DIMENSION IS 2 FOR SPIN POLARIZED
            dos_array = np.zeros((2,self.length))
            for m in range(6,self.length+6):
                line = dosstr[m].split()
                dos_array[0,count] = float(line[1])
                dos_array[1,count] = float(line[2])
                count+=1
        else:
            #FIRST DIMENSION IS 1 ELSE
            dos_array = np.zeros((1,self.length))
            for m in range(6,self.length+6):
                line = dosstr[m].split()
                dos_array[0,count] = float(line[1])
                count+=1

        return dos_array
    
    def energy(self):
        #READ ENERGIES FROM DOSCAR, RETURN AS ARRAY OF FLOATS
        dosstr = self.dosstr.split('\n')
        energy = np.zeros((1,self.length))
        count = 0
        for i in range(6,self.length+6):
            line = dosstr[i].split()
            energy[0,count] = float(line[0])
            count+=1
            
        return energy
    
    def odos(self):
        """ There's a lot here - the structure of DOSCAR for orbital case is quite complex. Take care 
            if modifying this section."""
        #HEADER INFO
        header = 7
        dosstr = self.dosstr.split('\n')
        start = header+self.length
        stop = len(dosstr)-1
        total_length = stop-start
        num_atoms = total_length/self.length
        energy = np.zeros((1,total_length))
        row_delete=0
        dividend=self.length
        
        #GET ENERGY
        count=0
        for n in range(start,stop):
            energy[0,count] = float(dosstr[n].split()[0])
            count+=1
        
        #CHECK FOR F_BANDS - VASP REQUIRES AN ADDITIONAL LINE IN SOC CASE
        if not energy[0,0]==energy[0,self.length+1]:
            fbands = True
            energy = energy[:,0:2*self.length]
            energy = energy[:,::2]
        else:
            fbands = False
            energy = energy[:,0:self.length]
        
        #GET DOS AS A FUNCTION OF CHARACTER FOR SPINPOL / NON_SOC
        if self.spinpol==True and self.soc==False:
            for test in range(start,start+1):
                line = dosstr[test].split()
                line = [float(j) for j in line]
            dos_array = np.zeros((2,total_length,(len(line)-1)/2))
            row_count=0
            for i in range(start,stop):                
                line = dosstr[i].split()
                if (i-start)%(dividend)==0 and i!=start:
                    row_delete+=1
                    dividend=(i-start)+self.length+1
                    continue
                line = [float(j) for j in line]
                col_count=0
                for k in range(1,len(line),2):
                    dos_array[0,row_count,col_count] = line[k]
                    dos_array[1,row_count,col_count] = line[k+1]
                    col_count+=1
                row_count+=1                
        
        #GET DOS AS A FUNCTION OF CHARACTER FOR NON_SPINPOL / NON_SOC
        if self.spinpol==False and self.soc==False:
            for test in range(start,start+1):
                line = dosstr[test].split()
                line = [float(j) for j in line]
            dos_array = np.zeros((1,total_length,(len(line)-1)))
            row_count=0
            for i in range(start,stop):
                if (i-start)%(dividend)==0 and i!=start:
                    row_delete+=1
                    dividend=(i-start)+self.length+1
                    continue                
                line = dosstr[i].split()
                line = [float(j) for j in line]
                col_count=0
                for k in range(1,len(line)):
                    dos_array[0,row_count,col_count] = line[k]
                    col_count+=1
                row_count+=1
                if row_count%self.length==0 and not row_count==0:
                    row_delete+=1
                    continue                
            
        #GET DOS AS A FUNCTION OF CHARACTER FOR SOC
        if self.soc==True and fbands==False:
            for test in range(start,start+1):
                line = dosstr[test].split()
                line = [float(j) for j in line]
            dos_array = np.zeros((4,total_length,(len(line)-1)/4))
            row_count=0
            for i in range(start,stop):
                if (i-start)%(dividend)==0 and i!=start:
                    row_delete+=1
                    dividend=(i-start)+self.length+1
                    continue                
                line = dosstr[i].split()
                line = [float(j) for j in line]
                col_count = 0
                for k in range(1,len(line),4):
                    dos_array[0,row_count,col_count] = line[k]
                    dos_array[1,row_count,col_count] = line[k+1]
                    dos_array[2,row_count,col_count] = line[k+2]
                    dos_array[3,row_count,col_count] = line[k+3]
                    col_count+=1
                row_count+=1
                if row_count%self.length==0 and not row_count==0:
                    row_delete+=1
                    continue                
        
        #GET DOS AS A FUNCTION OF CHARACTER WITH SOC AND FBANDS
        if self.soc==True and fbands==True:
            for test in range(start,start+1):
                line = dosstr[test].split()
                line = [float(j) for j in line]
            dos_array = np.zeros((4,total_length,(len(line)-1)/4+7))
            line_count=0
            row_count=0
            for i in range(start,stop):                
                if line_count%2==0:
                    if (i-start)%(dividend)==0 and i!=start:
                        print dosstr[i].split()
                        row_delete+=1
                        dividend=(i-start)+2*self.length+1
                        continue                    
                    line = dosstr[i].split()
                    line2 = dosstr[i+1].split()
                    total_line = np.append(line,line2,axis=0)
                    col_count=0
                    for k in range(1,len(total_line),4):
                        dos_array[0,row_count,col_count] = float(total_line[k])
                        dos_array[1,row_count,col_count] = float(total_line[k+1])
                        dos_array[2,row_count,col_count] = float(total_line[k+1])
                        dos_array[3,row_count,col_count] = float(total_line[k+2])
                        col_count+=1
                    row_count+=1
                    line_count+=1
                else:
                    line_count+=1
                if row_count%self.length==0 and not row_count==0:
                    row_delete+=1
                    continue                
                    
        #SUM FOR NUMBER OF ATOMS
        curr_length=dos_array.shape[1]
        new_length=curr_length-row_delete
        dos_array = dos_array[:,0:new_length]
        new_array = np.zeros((dos_array.shape[0],self.length,dos_array.shape[2]))
        for atom in range(num_atoms-1):
            new_array+=dos_array[:,self.length*atom:self.length*(atom+1),:]
        
        return energy,new_array             
                

#============================    DEFINE PROCAR CLASS   ==================================== 
class procar:
    """ Reader for VASP PROCAR file. This file is only generated for LORBIT tag being set.
        This file specifically provides orbital character information, and the file 
        structure can be complex. There is a more straighforward way to write these
        routines for the interested (I'm sure). But these should do the job."""
    
    def __init__(self,dir_procar,SpinPol=True,SOC=True):
        #READ IN PROCAR AS STRING
        self.dir_procar = dir_procar
        prostring = open(self.dir_procar+'PROCAR','r')
        prostring2 = prostring.read()
        prostring.close()
        self.prostring = prostring2
        #READ IN STRING OF ENERGIES
        energystring = os.popen('grep "# of bands:" '+self.dir_procar+'PROCAR')
        energystring2 = energystring.read()
        start = energystring2.find("# of bands:")
        stop = energystring2.find("# of ions:")
        bands = [int(s) for s in energystring2[start:stop].split() if s.isdigit()]
        self.nbands = bands[0]
        #READ IN KPOINT STRING
        kpointstring = os.popen('grep "k-point " '+self.dir_procar+'PROCAR').read()
        kpointstring2 = kpointstring.split('\n')
        kpointstring3 = filter(None,kpointstring2)  
        self.spinpol=SpinPol
        self.soc=SOC
        if self.spinpol==True and self.soc==False:
            self.nkpts = len(kpointstring3)/2
        else:
            self.nkpts = len(kpointstring3)
        return
    
    def nbands(self):
        #THIS IS CLEARLY NOT NEEDED - JUST RETURN SELF.NBANDS.
        energystring = os.popen('grep "# of bands:" '+self.dir_procar+'PROCAR')
        energystring2 = energystring.read()
        start = energystring2.find("# of bands:")
        stop = energystring2.find("# of ions:")
        bands = [int(s) for s in energystring2[start:stop].split() if s.isdigit()]
        numbands = bands[0]
        return numbands
    
    def nkpts(self):
        #SAME - JUST RETURN SELF.NKPTS
        kpointstring = os.popen('grep "k-point " '+self.dir_procar+'PROCAR').read()
        kpointstring2 = kpointstring.split('\n')
        kpointstring3 = filter(None,kpointstring2)
        nkpts = len(kpointstring3)
        return nkpts
    
    def kpoints(self):
        #RETURN ARRAY OF KPOINTS AS FLOAT.
        kpointstring = os.popen('grep k-point '+self.dir_procar+'PROCAR').read()
        kpointstring2 = kpointstring.split('\n')
        kvec = np.zeros((3,self.nkpts))
        col_count=0
        for i in range(1,len(kpointstring2)):
            kstr = kpointstring2[i]
            start = kstr.find(':')+1
            stop = kstr.find('weight')
            kstr_new = kstr[start:stop].replace('-',' -')
            try:
                if kstr_new == "":
                    continue
                else:
                    kstr_new2 = map(float,kstr_new.split())
                    kvec[:,col_count] = kstr_new2
                    col_count+=1
            except ValueError:
                break    
        return kvec
    
    def energies(self):
        counter1 = 0
        energyvec=[]
        os.system('echo "============= BAND ENERGIES =============="')
        #DEFINE DIMENSIONS OF ENERGY BANDS SPIN POL VS NON SPIN POL
        if self.spinpol==True and self.soc==False:
            enmat = np.zeros((2,self.nkpts,self.nbands))
        else:
            enmat = np.zeros((1,self.nkpts,self.nbands))
        while True:
            counter1+=1
            #ADJUST SPACING OF TEXT
            if counter1<10:
                spacestr = "  "
            elif counter1>=10 and counter1<100:
                spacestr = " "
            elif counter1>=100:
                spacestr = ""
            #GREP ENERGIES FROM PROCAR FILE
            os.system('echo "band '+spacestr+str(counter1)+'" '+self.dir_procar+'PROCAR')
            energystring =  os.popen('grep "band '+spacestr+str(counter1)+'" '+self.dir_procar+'PROCAR').read()
            energylist = energystring.split('\n')
            if not energystring:
                os.system(' echo "band '+str(counter1)+' NOT found" ')
                break
            else:
                #ITERATE THROUGH EACH LINE - SKIP KPOINTS READ ENERGIES
                for k in range(len(energylist)):
                    enstr = energylist[k]
                    start = enstr.find('energy')
                    stop = enstr.find('# occ')
                    try:
                        enlist = enstr[start:stop].split()
                        num = float(enlist[1])
                        if k==0:
                            energyvec = num
                        else:
                            energyvec = np.append(energyvec,num)
                    except IndexError:
                        continue
                if len(energyvec)==2*self.nkpts:
                    enmat[0,:,counter1-1] = energyvec[0:self.nkpts]
                    enmat[1,:,counter1-1] = energyvec[self.nkpts:]
                else:
                    enmat[0,:,counter1-1] = energyvec
        return enmat
    
    def labels(self):
        #READ LABELS - I.E. BAND CHARACTERS
        labelfile = os.popen('grep -m 2 "ion" '+self.dir_procar+'PROCAR').read()
        labelstring = labelfile.split('\n')[1]
        labels = map(str,labelstring.split())
        del labels[0]
        del labels[-1]
        return labels
    
    def character(self):
        nkpts = self.nkpts
        numbands = self.nbands
        soc = self.soc
        spinpol = self.spinpol
        start = [m.start() for m in re.finditer('ion ', self.prostring)]
        stop = [m.start() for m in re.finditer('band ', self.prostring)]
        nkcount=0
        nbcount=0
        os.system('echo "=============  BAND CHARACTERS  =============="')
        if self.spinpol==True and self.soc==False:
            for j in range(len(stop)):
                if j==(len(stop)-1):
                    segments = self.prostring[start[j]:].split("\n\n")
                else:
                    segments = self.prostring[start[j]:stop[j+1]].split("\n\n")
                count=0
                kpt=False
                for i in range(len(segments)):
                    if segments[i]:
                        try:
                            count+=1
                            block = segments[i].split("\n")
                            listblock = block[len(block)-1].split()
                            spinor = np.array([float(k) for k in listblock[1:len(listblock)-1]])
                            if count==1:
                                spinor1=spinor
                            #elif count==2:
                            #    spinor2=spinor
                        except ValueError:
                            os.system('echo "k-point '+str(nkcount+1)+'" '+self.dir_procar)
                            nkcount+=1
                            nbcount=0
                            kpt=True
                            break
                    else:
                        continue
                if j==0:
                    orbital_spin1 = np.zeros((nkpts,numbands,len(spinor1)))
                    orbital_spin2 = np.zeros((nkpts,numbands,len(spinor1)))
                    
                if kpt==False:
                    if nkcount==self.nkpts:
                        nkcount=0
                        nbcount=0
                    if nkcount<self.nkpts:
                        orbital_spin1[nkcount,nbcount,:] = spinor1
                        nbcount+=1
                    else:
                        print j, nkcount, nbcount
                        orbital_spin2[nkcount,nbcount,:] = spinor1
                        nbcount+=1
                
                elif kpt==True:
                    continue
            os.system('echo "k-point '+str(nkcount+1)+'" '+self.dir_procar)
            orbital = np.zeros((2,self.nkpts,self.nbands,orbital_spin1.shape[2]))
            orbital[0,:,:,:] = orbital_spin1
            orbital[1,:,:,:] = orbital_spin2 
        elif self.spinpol==False and self.soc==False:
            for j in range(len(stop)):
                if j==(len(stop)-1):
                    segments = self.prostring[start[j]:].split("\n\n")
                else:
                    segments = self.prostring[start[j]:stop[j+1]].split("\n\n")
                count=0
                kpt=False
                for i in range(len(segments)):
                    if segments[i]:
                        try:
                            count+=1
                            block = segments[i].split("\n")
                            listblock = block[len(block)-1].split()
                            spinor = np.array([float(k) for k in listblock[1:len(listblock)-1]])
                            if count==1:
                                spinor1=spinor
                            #elif count==2:
                            #    spinor2=spinor
                        except ValueError:
                            os.system('echo "k-point '+str(nkcount+1)+'" '+self.dir_procar)
                            nkcount+=1
                            nbcount=0
                            kpt=True
                            break
                    else:
                        continue
                if j==0:
                    orbital_spin1 = np.zeros((nkpts,numbands,len(spinor1)))
                    
                if kpt==False:
                    if nkcount==self.nkpts:
                        nkcount=0
                        nbcount=0
                    if nkcount<self.nkpts:
                        orbital_spin1[nkcount,nbcount,:] = spinor1
                        nbcount+=1
                    else:
                        print j, nkcount, nbcount
                        orbital_spin2[nkcount,nbcount,:] = spinor1
                        nbcount+=1
                
                elif kpt==True:
                    continue
            os.system('echo "k-point '+str(nkcount+1)+'" '+self.dir_procar)
            orbital = np.zeros((1,self.nkpts,self.nbands,orbital_spin1.shape[2]))
            orbital[0,:,:,:] = orbital_spin1             
        else:
            for j in range(len(stop)):
                if j==(len(stop)-1):
                    segments = self.prostring[start[j]:].split("\n\n")
                else:
                    segments = self.prostring[start[j]:stop[j+1]].split("\n\n")
                count=0
                kpt=False
                for i in range(len(segments)):
                    if segments[i]:
                        try:
                            count+=1
                            block = segments[i].split("\n")
                            listblock = block[len(block)-1].split()
                            spinor = np.array([float(k) for k in listblock[1:len(listblock)-1]])
                            if count==1:
                                spinor1=spinor
                            elif count==2:
                                spinor2=spinor
                            elif count==3:
                                spinor3=spinor
                            elif count==4:
                                spinor4=spinor
                        except ValueError:
                            os.system('echo "k-point '+str(nkcount+1)+'" '+self.dir_procar)
                            nkcount+=1
                            nbcount=0
                            kpt=True
                            break
                    else:
                        continue
                if j==0:
                    orbital_spin1 = np.zeros((nkpts,numbands,len(spinor1)))
                    orbital_spin2 = np.zeros((nkpts,numbands,len(spinor2)))
                    orbital_spin3 = np.zeros((nkpts,numbands,len(spinor3)))
                    orbital_spin4 = np.zeros((nkpts,numbands,len(spinor4)))
                    
                if kpt==False:
                    orbital_spin1[nkcount,nbcount,:] = spinor1
                    orbital_spin2[nkcount,nbcount,:] = spinor2
                    orbital_spin3[nkcount,nbcount,:] = spinor3
                    orbital_spin4[nkcount,nbcount,:] = spinor4
                    nbcount+=1
                
                elif kpt==True:
                    continue
            os.system('echo "k-point '+str(nkcount+1)+'" '+self.dir_procar)
            orbital = np.zeros((4,self.nkpts,self.nbands,orbital_spin1.shape[2]))
            orbital[0,:,:,:] = orbital_spin1
            orbital[1,:,:,:] = orbital_spin2
            orbital[2,:,:,:] = orbital_spin3
            orbital[3,:,:,:] = orbital_spin4
        
        return orbital
    
#===========================    DEFINE KPOINTS CLASS   ====================================
class kpoints:
    """ Read the kpoints file for information"""
    
    def __init__(self,dir_kpts):
        self.dir_kpts = dir_kpts
        kptfile = open(self.dir_kpts+'KPOINTS','r')
        kptstr = kptfile.read()
        kptfile.close()
        self.kptstr = kptstr
        return
        
    def kpath(self):
        #MEANT FOR SPECIAL POINTS FOR LINE MODE
        start = self.kptstr.find("rec")
        kfile = self.kptstr[start:].split("\n")
        path = len(kfile)
        kpath = np.zeros((3,(path-1)/3+1))
        k=0
        for j in range(1,path+1,3):
            try:
                kfile3 = kfile[j].split()
                kpath[0,k] = float(kfile3[0])  
                kpath[1,k] = float(kfile3[1])
                kpath[2,k] = float(kfile3[2])
                k= k+1        
            except IndexError:
                kfile3 = kfile[j-2].split()
                kpath[0,k] = float(kfile3[0])  
                kpath[1,k] = float(kfile3[1])
                kpath[2,k] = float(kfile3[2])
                break
        return kpath
    #ERROR IF NOT LINE MODE
    
    def nsample(self):
        #NEED TO COMPLETE
        pass
    
    def labels(self):
        start = self.kptstr.find("rec")
        lines = self.kptstr[start:].split('\n')
        kpoints=[]
        for i in range(1,len(lines)):
            try:
                sublines = lines[i].split()
                kpoints = np.append(kpoints,sublines[3])
            except IndexError:
                continue
        #Reprocess Kpoints
        for i in range(len(kpoints)):
            kstring = kpoints[i]
            newstring = kstring[1:2]
            kpoints[i] = newstring
        #Screen for duplicates
        klabel=[]
        for i in range(0,len(kpoints),2):
            try:
                if kpoints[i]=='G':
                    klabel = np.append(klabel,'$\Gamma$')
                else:
                    klabel = np.append(klabel,kpoints[i])
            except IndexError:
                break
        klabel = np.append(klabel,kpoints[len(kpoints)-1])
        return klabel 

class qscript:
    """ Read information from marcc.job qscript file"""
    def __init__(self,filepath):
        self.dir_qscript = filepath
        qstr = open(self.dir_qscript+'marcc.job','r')
        qstr2 = qstr.read()
        qstr.close()
        self.qstr = qstr2
        
    def partition(self):
        start = self.qstr.find('partition')
        stop = self.qstr.find('module')
        queuestr2 = self.qstr[start:stop].split('\n')[0]
        start2 = queuestr2.find('=')
        partition = queuestr2[(start2+1):]
        return partition
    
    def time(self):
        start = self.qstr.find('time')
        stop = self.qstr.find(':0:0')
        queuestr2 = self.qstr[start:stop].split('\n')[0]
        start2 = queuestr2.find('=')
        time = queuestr2[(start2+1):]
        time = time*3600
        return time
        
class findreplace:
    def __init__(self,filename,filepath,searchphrase,replacephrase):
        self._filepath = filepath
        self._filename = filename
        self._searchphrase = searchphrase
        self._replacephrase = replacephrase
        oldfile=open(self._filepath+self._filename,'r')
        oldstring=oldfile.read()
        oldfile.close()   
        newstring=oldstring.replace(searchphrase,replacephrase)   
        newfile=open(self._filepath+self._filename,'w')
        newfile.write(newstring)
        newfile.close()
        return
    
#===========================    DEFINE INCAR CLASS   =====================================
class incar:
    """ Read information from INCAR file - specifically soc and spinpol setting."""
    def __init__(self,dir_incar):
        self.dir_incar = dir_incar
        
    def soc(self):
        soc_string = os.popen("grep LSORBIT "+self.dir_incar+"INCAR")
        soc_string = soc_string.read()
        start = soc_string.find('=')
        stop = soc_string.find('\n')
        start_new = start+1
        stop_new = stop-1
        if soc_string[start_new:stop_new]==".FALSE." or soc_string[start_new:stop_new]==" .FALSE." or soc_string[start_new:stop_new]==" .FALSE. ":
            soc = bool(False)
        elif soc_string[start_new:stop_new]==".TRUE." or soc_string[start_new:stop_new]==" .TRUE." or soc_string[start_new:stop_new]==" .TRUE. ":
            soc = bool(True)
        else:
            soc = bool(False) 
            
        return soc
    
    def spin(self):
        spin_string = os.popen("grep ISPIN "+self.dir_incar+"INCAR")
        spin_string = spin_string.read()
        start = spin_string.find('=')+1
        stop = spin_string.find('\n')
        if spin_string[start:stop]=="1" or spin_string[start:stop]==" 1" or spin_string[start:stop]==" 1 ":
            spin = bool(False)
        elif spin_string[start:stop]=="2" or spin_string[start:stop]==" 2" or spin_string[start:stop]==" 2 ":
            spin = bool(True)
        else:
            spin = bool(False) 
            
        return spin
    

#==================================            ===========================================
#==================================    END     ===========================================
#==================================            ===========================================