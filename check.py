## JAKE A TUTMAHER
## JOHNS HOPKINS UNIVERSITY
## MCQUEEN LABORATORY
## THE INSTITUTE FOR QUANTUM MATTER
## DEPARTMENT OF PHYSICS, DEPARTMENT OF CHEMISTRY, DEPARTMENT OF MATERIALS SCIENCE AND ENGINEERING
##
## CONTACT: jtutmah1@jhu.edu
###################################################################################################

""" Iterate through subdirectories which contain various VASP jobs and return information such as
    convergence and computational statistics.
"""
#==============================  IMPORT MODULES  ===============================================
import os 
import glob
import numpy as np
from collections import OrderedDict
import csv
import IQM.VASPread as VASPread
import utils.inputs as inputs

#=============================      GET SUBDIRS     ============================================
#GET CURRENT WORKING DIRECTORY
maindir = inputs.get_current_directory()

#GET ALL SUBDIRECTORIES FOR VASP RUNS
names = inputs.get_immediate_subdirectories(maindir)
length = len(names)

#=============================   DICT HEADINGS   ===============================================
#GENERATE DICT KEYS
columnName0='COMPOUND'
columnName1='RELAX'
columnName2='NUM RUNS'
columnName3='RELAX TIME (HRS)'
columnName4='REX NUM VALENCE'
columnName5='STATIC'
columnName6='STATIC TIME (HRS)'
columnName7='STAT NUM VALENCE'
columnName8='DFPT'
columnName9='DFPT TIME (HRS)'
columnName10='DFPT NUM VALENCE'

d = OrderedDict()

#=================================  KEYPHRASES  =================================================
#THESE PHRASE ARE AT THE END OF EACH STD OUT (SLURM.####) FILE IF EACH RUN RAN TO COMPLETION
strallow=' reached required accuracy - stopping structural energy minimisation'
strallow2=' writing wavefunctions'
strallow3=' Linear response finished'

#==============================  CHECK FUNCTION  ================================================
def check(dir_main,prefix,keyphrase,subdir_length,subdir):
    """ This function iterates through the std out files in each directory, finds the most
        recent one, and checks if it's been run to completion or not. @dir_main is the main 
        directory containing these jobs, @prefix is the specific run type (i.e. relax, static,
        dfpt) and keyphrase is the std out keyphrase to find."""
    
    #CLEAR ARRAYS    
    Entries=[]
    Number = [] 
    
    #ITERATE THROUGH ALL SUBDIRECTORIES
    for j in range(subdir_length):
        
        #DEFINE PATH OF SPECIFIC DIRECTORY
        dir_out = dir_main+subdir[j]+'/'+prefix+'/'
        path = dir_out
        files = dir_out+'*.out'
        
        #IF THIS PATH ACTUALLY CONTAINS A RELAX, STATIC, OR DFPT DIR
        if os.path.exists(dir_out):
            
            #LIST ALL .OUT FILES AS ARRAY
            name = glob.glob(files)
            name = np.asanyarray(name)
            
            #IF THERE ARE SLURM FILES, LOOP THROUGH AND SEARCH FOR KEYWORD
            if len(name)!=0:
                Number = np.append(Number,len(name))
                num_vec = []
                
                #PULL ID NUMBER FOR ALL .OUT FILES CONTAINED IN DIRECTORY
                for k in range(len(name)):
                    name2 = name[k]
                    num = int(name2[-11:-4])
                    num_vec = np.append(num_vec,num)
                    
                #FIND .OUT FILE WITH MAX NUMBER (MOST RECENT NUMBER) AND READ AS STRING
                m = max(num_vec)
                position = [i for i, j in enumerate(num_vec) if j == m]
                str_output = os.popen('grep "'+ keyphrase +'" '+name[position][0])
                string = str_output.read()
                
                #IF KEYPHRASE EXISTS FROM GREP - THEN IT HAS CONVERGED
                if string:
                    Entries=np.append(Entries,' Y ')
                else:
                    Entries=np.append(Entries,' N ')
        #OUTPUT FILES NOT FOUND            
            else:
                Entries=np.append(Entries,' DNR ')
                Number=np.append(Number,0)
        else:
            Entries=np.append(Entries,'DNR')
            Number=np.append(Number,0)
            
    return Entries,Number

#====================================  CHECK OUTPUT  ==========================================

[Entries1, Number] = check(maindir,'RELAX',strallow,length,names)
[Entries2,Number2] = check(maindir,'STATIC',strallow2,length,names)
[Entries3,Number3] = check(maindir,'DFPT',strallow3,length,names)

#=====================================  READ TIME  =============================================
timerex = []
timestat = []
timedfpt = []

for k in range(length):
    
    #DEFINE RELAX, STATIC, AND DFPT DIRECTORIES
    relax_dir = maindir+names[k]+'/RELAX/'
    static_dir = maindir+names[k]+'/STATIC/'
    dfpt_dir = maindir+names[k]+'/DFPT/'
    
    #TRY READING TIME - THROW EXCEPTION IF FILE DOES NOT EXIST
    try:
        relax_outcar = VASPread.outcar(relax_dir)
        rex_time = relax_outcar.time()
    except (IOError,IndexError):
        rex_time = 'DNR'
    try:
        static_outcar = VASPread.outcar(static_dir)
        stat_time = static_outcar.time()
    except (IOError,IndexError):
        stat_time = 'DNR'
    try:
        dfpt_outcar = VASPread.outcar(dfpt_dir)
        dfpt_time = dfpt_outcar.time()
    except (IOError,IndexError):
        dfpt_time = 'DNR'
    
    #APPEND RELAX, STATIC, AND DFPT TIME TO ARRAY
    timerex = np.append(timerex,rex_time)
    timestat = np.append(timestat,stat_time)
    timedfpt = np.append(timedfpt,dfpt_time)

#=================================  READ VALENCE  ===========================================
rex_valence = []
stat_valence = []
dfpt_valence = []

for l in range(length):
    #DEFINE RELAX, STATIC, AND DFPT DIRECTORIES
    relax_dir = maindir+names[l]+'/RELAX/'
    static_dir = maindir+names[l]+'/STATIC/'    
    dfpt_dir = dfpt_dir = maindir+names[l]+'/DFPT/'
    
    #PULL NUMVALENCE FROM EIGENVAL - EXCEPT IF FILES ARE EMPTY OR DON'T EXIST
    try:
        rex_eigenval = VASPread.eigenval(relax_dir)
        rex_nvalence = rex_eigenval.nvalence()
    except (IOError,IndexError):
        rex_nvalence = 'DNR'       
    try:
        stat_eigenval = VASPread.eigenval(static_dir)
        stat_nvalence = stat_eigenval.nvalence()
    except (IOError,IndexError):
        stat_nvalence = 'DNR'    
    try:
        dfpt_eigenval = VASPread.eigenval(dfpt_dir)
        dfpt_nvalence = dfpt_eigenval.nvalence()
    except (IOError,IndexError):
        dfpt_nvalence = 'DNR'   
    
    #WRITE NVALENCE TO ARRAY
    rex_valence = np.append(rex_valence,rex_nvalence)
    stat_valence = np.append(stat_valence,stat_nvalence)
    dfpt_valence = np.append(dfpt_valence,dfpt_nvalence)
    
#======================   READ AND PRINT SUMMARY   ================================
#DEFINE DICTIONARY ENTRIES
d[columnName0]=names
d[columnName1]=Entries1
d[columnName2]=Number
d[columnName3]=timerex
d[columnName4]=rex_valence
d[columnName5]=Entries2
d[columnName6]=timestat
d[columnName7]=stat_valence
d[columnName8]=Entries3
d[columnName9]=timedfpt
d[columnName10]=dfpt_valence

#WRITE HEADER INFORMATION
header = d.keys()
writer = csv.writer(open(maindir+'summary.csv', 'wb'))
writer.writerow(header)
writer.writerows(zip(*d.values()))
