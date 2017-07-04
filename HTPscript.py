
#===============================  MODULES  ===================================
import os
import IQM.manage as manage

#=============================  DIRECTORIES  =================================
dir1 = os.popen("pwd").readlines()[0]
length2=len(dir1)
dir1 = dir1[:length2-1]+'/'

#=======================  DEFINE HTP ENVIRONMENT  =============================
dir_main = dir1
dir_cifs = dir1+'[PATH TO CIFS]'
dir_temp = '[PATH TO TEMPLATE FILES FOR VASP]'
dir_pot = '[PATH TO VASP POTCARS]'

sampleSet1 = manage.HTP(dir_main,dir_cifs,dir_temp,dir_pot,'sbatch marcc.job')

#===========================  CONVERGE RELAX =================================
print('BEGINNING RELAXATION')
sampleSet1.convergeAll('RELAX','RELAX',1000)

#=======================  RUN STATIC FROM RELAX ===========================
print('BEGINNING STATIC')
sampeSet1.continueAll('RELAX','STATIC','STATIC')
print('HTP SCRIPT TERMINATED')






        
        

        
            
  
