## JAKE A TUTMAHER
## JOHNS HOPKINS UNIVERSITY
## MCQUEEN LABORATORY
## THE INSTITUTE FOR QUANTUM MATTER
## DEPARTMENT OF PHYSICS, DEPARTMENT OF CHEMISTRY, DEPARTMENT OF MATERIALS SCIENCE AND ENGINEERING
##
## CONTACT: jacob.tutmaher@gmail.com
###################################################################################################
""" Heigh Throughput management algorithms for MARCC, specifically designed for VASP.
    These classes generate, monitor, and converge jobs automatically by calling 
    HTP in terminal.
    
    These classes are implemented in PYTHON/HTPscript.py
    
    ---
    
"""

#===========================    IMPORT SPECIFIC PACKAGES   ============================
import os
import glob
import numpy as np
import VASPread
import time

#===========================    DEFINE HTP CLASS   ====================================
class HTP: 
    """ HTP (High Throughput) is a class that uses a given set of cifs (@self.cif_dir) and 
        TEMPLATES (@self.temp_dir) - along with the Potential File directory (@self.pot_dir)
        to generate an arbitrary number of jobs (equal to the number of cifs in the cif director)
        and complete them to convergence. 
        
        Note: This class is implemented in PYTHON/HTPscript.py
    """
    def __init__(self,dir1,dir2,dir3,dir4,systemcommand):
        """ Initialize local variables : @self.main_dir (parent directory from which all jobs 
            are executed), @self.cif_dir (location of all the cif files), and @self.temp_dir
            and @self.pot_dir which should be set on MARCC. The system command is sbatch 
            marcc.job, and the keyphrase which is the last line of the std out file.
        """
        self.main_dir = dir1
        self.cif_dir = dir2        
        self.temp_dir = dir3
        self.pot_dir = dir4
        self.systemcommand = systemcommand
        self.keyphrase = ' reached required accuracy - stopping structural energy minimisation'
        
    def cifs(self):
        """ Return a list of the cif names located in @self.cif_dir.
            Note, cif files are typically crystal structure files
        """
        os.chdir(self.cif_dir)
        names = glob.glob(os.path.basename('*.cif'))
        for i in range(len(names)):
            names[i]=os.path.splitext(names[i])[0]
        self._cifs = names
        return self._cifs
        
    def number(self):
        """ Return the number of cifs in @self.cif_dir
        """
        self._number = len(self.cifs())
        return self._number
    
    def runAll(self,tempprefix,finalprefix):
        """ Run all compounds in a given folder by generating jobs from the @self.cif_dir cifs.
            Uses the STP class given below.
        """
        system = STP(self.main_dir,self.cif_dir,self.temp_dir,self.pot_dir,self.systemcommand)
        #GENERATE JOBS
        for i in range(self.number()):
            system.genJob(tempprefix,finalprefix,self._cifs[i])
        #SUBMIT JOBS
        for j in range(self.number()):
            status = system.runJob(self._cifs[j],finalprefix,self.keyphrase)
        return
    
    def logfile(self,jobtype,numrun,numconverge,numresub,iteration):
        """ Create a logfile (log.txt) for a given job type. Implemented in methods listed below.
        """
        logfile=open(self.main_dir+'log.txt','a')
        titlestr = '=========================  '+jobtype+' ITERATION '+str("{0:.2f}".format(iteration))+' HRS  ======================='
        logfile.write(titlestr+'\n'+' Number Submitted: '+str(numresub)+'\n Number Running: '+str(numrun)+'\n Number Converged: '+str(numconverge)+'\n \n')
        logfile.close()
        return
            
    def convergeAll(self,tempprefix,finalprefix,interval):
        """ This method is used for relaxation runs - which may take a few iterations to converge. 
            It generates the initial relaxation jobs and then runs them until convergence. 
            @tempprefix: the subdirectory name for the TEMPLATE INCAR file.
            @finalprefix: the subdirectory you will copy the INCAR file to.
            This method runs for a maximum time of tmax, and exits the while loop when all 
            jobs have been converged.
            
            This method is implemented in PYTHON/HTPscripty.py for Relaxation jobs.
        """
        #INITIALIZE STANDARD THROUGHPUT FOR EACH JOB
        system = STP(self.main_dir,self.cif_dir,self.temp_dir,self.pot_dir,self.systemcommand)
        t=0
        counter = 1
        
        #READ QFILE FOR PARTITION INFORMATION
        qfile = VASPread.qscript(self.temp_dir+tempprefix+'/')
        partition = qfile.partition()
        
        #MAY NEED ADJUSTED DEPENDING ON TIME CONSTRAINTS
        tmax = 259000
        
        #ENTER WHILE LOOP FOR CONVERGENCE
        while t<=tmax:
            #SET TO ZERO - AS JOBS WERE JUST SUBMITTED
            numrun = 0
            numconverge = 0
            numresub = 0
            
            #READ LIST OF JOBS FOR THAT PARTITION
            output = os.popen('squeue -p '+partition+' -u jtutmah1@jhu.edu').readlines()
            name=None
            
            #ENTER EACH SUBDIRECTORY THAT IS RUNNING JOBS
            for i in range(self.number()):
                finaldir = self.main_dir+self._cifs[i]+'/'+finalprefix+'/'
                if os.path.exists(finaldir):
                    run = 1
                    #SEE IF NAME OF COMPOUND IS STILL IN QUEUE
                    for j in range(1,len(output)):
                        name = output[j].split()[2]
                        try:    
                            stop = name.index('_')
                            name = name[0:stop]
                            cif_name = self._cifs[i]
                        except ValueError:
                            stop = len(name)                
                            name = name[0:stop]
                            cif_name = self._cifs[i] #CHANGE CIF NAME LENGTH TO MATCH FIRST CHARACTERS
                            cif_name = cif_name[0:stop]
                        if cif_name==name:
                            run=0
                            break
                        else:
                            run=1
                    #IF NAME NOT IN QUEUE - RESUBMIT, UNLESS CONVERGED, THEN CONTINUE
                    if run==1:
                        converge = system.runJob(self._cifs[i],finalprefix,self.keyphrase)
                        if converge == True:
                            numconverge += 1
                        else:
                            os.system('echo resubmitting')                            
                            numresub += 1
                        
                    else:
                        numrun += 1
                        continue
                #IF DIRECTORY WAS FOUND TO NOT EXIST, THEN GENERATE AND SUBMIT JOB (FIRST PASS)
                else:
                    system.genJob(tempprefix,finalprefix,self._cifs[i])
                    time.sleep(2)
                    string = system.runJob(self._cifs[i],finalprefix,self.keyphrase)
                    numresub += 1
            
            #TIME IN MINUTES FOR LOG FILE
            iteration = t/3600.0
            self.logfile(finalprefix,numrun,numconverge,numresub,iteration)
            
            #IF ALL CIFS ARE CONVERGED
            if numconverge==self.number() and not counter==1:
                t = tmax+1
                
            #ELSE CONTINUE IN WHILE LOOP
            else:
                counter = counter+1
                t = t+interval
                time.sleep(interval)        
        os.system('echo CONVERGE FINISHED')
        return
        
    def dependentAll(self,tempprefix,finalprefix,dependentdir):
        """ This method is currently unused - and incomplete. It's meant to provide dependent checks. 
            A similar method is currently implemeneted in continueAll.
        """
        system=STP(self.main_dir,self.cif_dir,self.temp_dir,self.pot_dir,self.systemcommand)
        for i in range(self.number()):
            tempdir = self.temp_dir+tempprefix
            finaldir = self.main_dir + self._cifs[i]+finalprefix
            ##Check Log File
            
    def continueAll(self,initprefix,finalprefix,tempprefix):
        """ Continue all sets up a continuation run for all subdirectories in a main folder.
            This means that a new directory is created (@finalprefix), and the INCAR, POTCAR,
            and CONTCAR -> POSCAR are copied over from @initprefix. Whether it's a true 
            continuation or a restart depends on the settings in the INCAR (@tempprefix) 
            settings. See VASP documentation for information on how to implement a continuation run.
            
            This method is implemented in PYTHON/HTPscript.py for Static jobs.
        """
        #SET UP STP ENVIRONMENT
        system = STP(self.main_dir,self.cif_dir,self.temp_dir,self.pot_dir,self.systemcommand)
        numresub = 0
        
        #GENERATE CONTINUATION JOBS FOR ALL SUBDIRECTORIES - I.E. NUM SUBDIRECTORIES EQUALS NUM CIFS
        for i in range(self.number()):
            system.genContinueJob(initprefix,finalprefix,tempprefix,self._cifs[i])
            
        #RUN CONTINUATION JOBS FOR ALL SUBDIRECTORIES
        for j in range(self.number()):
            system.runJob(self._cifs[j],finalprefix,self.keyphrase)
            numresub += 1
        
        #PRINT INFO TO LOG FILE
        self.logfile(finalprefix,0,0,numresub,0)
        return
    
    def continueSuperAll(self,initprefix,finalprefix,tempprefix):
        """ This method is similar to continueAll - except that it specifically generates
            supercells from the original CONTCAR (from @initprefix). It sets up a STP 
            environment first, and then runs the genSuperJob. Please see the genSuperJob method
            in the STP class for more information on the role of valence electrons on 
            supercell size.
            
            This method is implemented in PYTHON/HTPscript.py for DFPT jobs.
        """
        #SET UP STP ENVIRONMENT 
        system = STP(self.main_dir,self.cif_dir,self.temp_dir,self.pot_dir,self.systemcommand)
        numresub = 0
        
        #ITERATE THROUGH ALL SUBDIRECTORIES - I.E. NUMBER OF CIFS
        for i in range(self.number()):
            initdir = self.main_dir+'/'+self._cifs[i]+'/'+initprefix+'/'
            #GET VALENCE ELECTRONS FROM NUM VALENCE
            try:
                eigenfile = VASPread.eigenval(initdir)
                try:
                    numvalence = eigenfile.nvalence()
                except IndexError:
                    numvalence = 1000
            except IOError:
                numvalence = 1000
            #GENERATE CONTINUATION JOB AND THEN SUPERCELL 
            system.genContinueJob(initprefix,finalprefix,tempprefix,self._cifs[i])
            system.genSuperJob(finalprefix,self._cifs[i],numvalence)
            
        #SUBMIT/RUN ALL JOBS
        for j in range(self.number()):
            system.runJob(self._cifs[j],finalprefix,self.keyphrase)
            numresub += 1
            
        #UPDATE LOG FILE
        self.logfile(finalprefix,0,0,numresub,0)
        return    
    
        
        
class STP:
    """ The STP class essentially implements the same methods as HTP - just at a single
        job level. This methods are commonly called iteratively in the HTP class above.
    """
    def __init__(self,dir1,dir2,dir3,dir4,systemcommand):
        """ Initialize variables: the parent directory (self.main_dir), the cif directory
            (@self.cif_dir), the template directory (@self.temp_dir), and the potential
            directory (@self.pot_dir). The system command is the command needed to 
            submit jobs on your cluster. On marcc, it is sbatch + qscript_name
        """
        self.main_dir = dir1
        self.cif_dir = dir2        
        self.temp_dir = dir3
        self.pot_dir = dir4
        self.systemcommand = systemcommand
        
    def genJob(self,tempprefix,finalprefix,cifname):
        """ Generates the five main scripts - INCAR, POTCAR, POSCAR, KPOINTS, and marcc.job
            needed to run a job in VASP. This scripts either pull information from a template
            directory (INCAR - @tempprefic), or from the cif (POSCAR and POTCAR - @cifname)
            It also uses a java kpoint generator developed by the Mueller Group to generate
            a kpoint file based on each POSCAR.
            
            This method is typically implemented for Relaxation runs - see convergeAll.
        """
        #DEFINE RELEVANT DIRECTORIES
        finaldir = self.main_dir+cifname+'/'+finalprefix+'/'
        tempdir = self.temp_dir+tempprefix+'/'
        if not os.path.exists(finaldir):
            os.makedirs(finaldir)
        
        #GENERATE INPUTS IN THE FINAL DIRECTORY
        Inputs = inputs(cifname)
        Inputs.INCAR(tempdir,finaldir,'INCAR')
        Inputs.POSCAR_POTCAR(self.cif_dir,finaldir)
        Inputs.KPOINTS(self.main_dir,finaldir)
        Inputs.QUEUE(tempdir,finaldir,finalprefix,'marcc.job')
        
        #MOVE BACK TO MAIN DIRECTORY
        os.chdir(self.main_dir)
        return

    def runJob(self,cifname,runprefix,keyphrase):
        """ This method runs the sbatch command to submit the given job (@cifname+@runprefix) 
            on the cluster. It enters the directory, and searches it to ensure that the job
            has not already converged. It finds the slurm.#### file with the highest number,
            and searches it for the converged @keyphrase. 
            
            This method is used in the HTP.convergeAll method above, along with other methods.
        """
        #SEARCH STD OUT FILES
        rundir = self.main_dir+cifname+'/'+runprefix+'/'
        path = rundir+'*.out' 
        name = glob.glob(path)
        name = np.array(name)
        
        #ITERATE THROUGH PREVIOUSLY RUN JOBS
        if not len(name) == 0 and not keyphrase==False:            
            num_vec = []
            for k in range(len(name)):
                name2 = name[k]
                num = int(name2[-11:-4])
                num_vec = np.append(num_vec,num)
                
            #FIND MOST RECENT STD OUT FILE
            m = max(num_vec)
            position = [i for i, j in enumerate(num_vec) if j == m]
            input_STDOUT = open(name[position],'r').read()
            STDOUT_str = input_STDOUT.split('\n')
            str_comp1 = STDOUT_str[-3]
            str_comp2 = STDOUT_str[-2]
            str_comp3 = STDOUT_str[-4]
            
            #SUBMIT ACCORDINGLY
            if str_comp1 == keyphrase or str_comp2 == keyphrase or str_comp3 == keyphrase:
                converge = True
            else:
                converge = False
                os.chdir(rundir)
                contcarfile = open('CONTCAR','r')
                contcarstring = contcarfile.read()
                if not contcarstring:
                    os.system(self.systemcommand)
                    os.chdir(self.main_dir)
                else:
                    os.system('cp CONTCAR POSCAR')
                    os.system(self.systemcommand)
                    os.chdir(self.main_dir)                
        
        #IF NO STD OUTS EXIST, THEN RUN THE JOB (FIRST PASS)
        else:
            converge = False
            os.chdir(rundir)
            os.system(self.systemcommand)
            os.chdir(self.main_dir)
        return converge
    
    def genContinueJob(self,initprefix,finalprefix,tempprefix,cifname):
        """ This generates a continuation job using the CONTCAR contained in @initprefix and
            the template files contained in @tempprefix. Note, that unlike STP.genJob, this 
            method AUTOMATICALLY implements a gamma-centered 9x9x9 kpt mesh. May be adjusted
            to implement an arbitrary k-mesh.
            
            This method is implemented in the HTP.continueAll method above.
        """
        #DEFINE AND CREATE RELEVANT DIRECTORIES
        initdir = self.main_dir+cifname+'/'+initprefix+'/'
        finaldir = self.main_dir+cifname+'/'+finalprefix+'/'
        tempdir = self.temp_dir+tempprefix+'/'
        if not os.path.exists(finaldir):
            os.makedirs(finaldir)
            
        #COPY OVER FILES (I.E. CONTCAR AND POTCAR) FROM OLD JOB TO NEW DIRECTORY
        os.system('cp '+initdir+'CONTCAR '+finaldir+'POSCAR')
        os.system('cp '+initdir+'POSCAR_original '+finaldir+'POSCAR_original')
        os.system('cp '+initdir+'POTCAR '+finaldir+'POTCAR')
        
        #INITIALIZE AND GENERATE INPUT FILES
        Inputs = inputs(cifname)
        Inputs.INCAR(tempdir,finaldir,'INCAR')
        
        #9X9X9 MESH IF STATIC JOB
        if finalprefix=="STATIC":
            os.chdir(finaldir)
            kpoints=open("KPOINTS",'w')
            kpoints.write('k-points\n 0\nGamma\n 9 9 9\n 0 0 0\n')
            kpoints.close()
            os.chdir(self.main_dir)
        else:
            Inputs.KPOINTS(self.main_dir,finaldir)
            
        #WRITE QUEUE SCRIPT
        Inputs.QUEUE(tempdir,finaldir,finalprefix,'marcc.job')
        os.chdir(self.main_dir)
        return
            
    def genSuperJob(self,finalprefix,cifname,numvalence):
        """ This generates a supercell (continuation) job using PHONOPY. Any unitcell with a valence count below
            75 will form a 2x2x2 supercell - otherwise it will form a 1x1x1 supercell. This limit can 
            of course be adjusted, this is just based on my experience with computation time. NOTE, this 
            method also uses a 9x9x9 kpoint mesh - rather than a kpoint generator. The kpt density here
            may also be adjusted in the code below as needed.
        """
        #FINAL - SUPERCELL DIRECTORY
        super_dir = self.main_dir+cifname+'/'+finalprefix
        
        #THE SUPERCELL DIRECTORY MUST BE CREATED FIRST USING GENCONTINUE - OTHERWISE THE POSCAR DOES NOT EXIST YET
        if not os.path.exists(super_dir):
            return 
        
        #POSCAR AVAILABLE - MAKE SUPERCELL JOB. NOTE, THE COMMANDS FOLLOW THOSE OUTLINED IN PHONOPY DOCUMENTATION
        else:
            os.chdir(super_dir)
            os.system('mv POSCAR POSCAR-unitcell')
            if (8*numvalence < 600):
                os.system('phonopy -d --dim="2 2 2" -c POSCAR-unitcell')
                os.system('mv SPOSCAR POSCAR')
                mesh = open('mesh.conf','w')
                mesh.write('ATOM_NAME = '+cifname+'\nDIM = 2 2 2\nDOS=.TRUE.\nTETRAHEDRA=.TRUE.\nMP = 31 31 31\nFORCE_CONSTANTS = READ\n')
                mesh.close
                band = open('band.conf','w')
                band.write('ATOM_NAME = '+cifname+'\nDIM = 2 2 2\nQPOINTS=.TRUE.\nFORCE_CONSTANTS=READ\n')
                band.close()
            else:
                os.system('phonopy -d --dim="1 1 1" -c POSCAR-unitcell')
                os.system('mv SPOSCAR POSCAR')
                mesh = open('mesh.conf','w')
                mesh.write('ATOM_NAME = '+cifname+'\nDIM = 1 1 1\nDOS=.TRUE.\nTETRAHEDRA=.TRUE.\nMP = 31 31 31\nFORCE_CONSTANTS = READ\n')
                mesh.close
                band = open('band.conf','w')
                band.write('ATOM_NAME = '+cifname+'\nDIM = 1 1 1\nQPOINTS=.TRUE.\nFORCE_CONSTANTS=READ\n')
                band.close() 
                
            #INITALIZE AND GENERATE NEW KPOINT FILE
            Inputs = inputs(cifname)
            #Inputs.KPOINTS(self.temp_dir,super_dir)
            #Modified Kpoints generation with following 3 lines
            kpoints=open("KPOINTS",'w')
            kpoints.write('k-points\n 0\nGamma\n 9 9 9\n 0 0 0\n')
            kpoints.close() 
            
            #RETURN TO MAIN DIRECTORY - DON'T LET MAIN SCRIPT STAY IN SUBDIR
            os.chdir(self.main_dir)
            return               
                
class inputs:
    """ This class contains methods to generate the 5 main input scripts used by VASP - INCAR,
        POSCAR, POTCAR, KPOINTS, and marcc.job (qscript). A few of the methods rely on java 
        routines contained in the JAVA folder. See java src file for more details.
    """
    def __init__(self,name):
        self.name = name
        
    def INCAR(self,tempdir,finaldir,incarname):
        #FIND TEMPLATE INCAR FILE, COPY OVER TO FINAL DIR AND REPLACE COMPOUND NAME
        os.system('cp '+tempdir+incarname+' '+finaldir+'INCAR')
        VASPread.findreplace('INCAR',finaldir,'$NAME$',self.name)
        return
    
    def POSCAR_POTCAR(self,cifpath,finaldir):
        #GENERATE POSCAR AND POTCAR FROM CIF USING THE GENPOSCAR JAVA PROGRAM (CREDIT PANDU)
        os.chdir(cifpath)
        os.system('cp '+self.name+'.cif '+finaldir)
        os.chdir(finaldir)
        os.system('genPoscar')
        os.system('mv POSCAR.vasp POSCAR')
        os.system('cp POSCAR POSCAR_original')
        #PHONOPY IS USED TO REGENERATE POSCAR IN SETTING CONSISTENT WITH PHONOPY SETTING
        os.system('phonopy --symmetry')
        os.system('cp PPOSCAR POSCAR')
        os.system('mv POTCAR.vasp POTCAR')
        return        
        
    def POSCAR(self,cifpath,finaldir):
        """ THIS METHOD IS OUTDATED - use POSCAR_POTCAR to generate these files."""
        os.chdir(cifpath)
        os.system('cif2cell '+self.name+'.cif -p vasp')
        if os.path.exists(cifpath+'POSCAR'):
            os.system('mv POSCAR '+finaldir)
            os.system('cp '+finaldir+'POSCAR '+finaldir+'POSCAR_original')
            return
        else:
            return  
        
    def POTCAR(self,potdir,finaldir):
        """ Generates POTCAR based on atoms in POSCAR. OUTDATED now."""
        if os.path.exists(finaldir+'POSCAR_original'):
            posfile = open(finaldir+'POSCAR_original','r')
            postring=posfile.read()
            newstring = postring.split('\n')
            newstring1 = newstring[0]
            start = newstring1.find('Species order:')
            newstring2 = newstring1[start:].split()
            Atoms=newstring2[2:]
            open(finaldir+'POTCAR','a').close()
            with open(finaldir+'POTCAR', 'w') as outfile:
                for fname in Atoms:
                    with open(potdir+fname) as infile:
                        outfile.write(infile.read())
            return
        else:
            return
                    
    def KPOINTS(self,tempdir,finaldir):
        """ Uses kpoint generator in GENKPTS folder."""
        precalc='/home-2/jtutmah1@jhu.edu/work/jake/GENKPTS/'
        os.system('cp '+precalc+'PRECALC '+finaldir)
        os.chdir(finaldir)
        os.system('genKPTS')
        os.chdir(tempdir)
        return
    
    def QUEUE(self,tempdir,finaldir,finalprefix,script):
        """ This method copies marcc.job script from template files and replaces the $NAME$ keyword
            with the actual compound name.
        """
        os.system('cp '+tempdir+script+' '+finaldir)
        VASPread.findreplace(script,finaldir,'$NAME$',self.name+'_'+finalprefix)
        return
    
    
    
        
    

        
    