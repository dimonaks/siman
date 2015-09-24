#!/usr/bin/env python
"""  
1. need anfuncs module
2. no argument-read from database lofl.s, 1 or more - read from real files

3. Rules of naming
first digit - type of calculation
 0 - cell with grain boundary 
 1 - cell without grain boundary
 2 - cell with carbon in position 605 (situated in dir puregb) (it is one exclusion from rules-must be deprecated)
600,700,800 - octahedral site without gb
 601-606
 701-706
 801-806 #see below
second digit - number of set.
"""
import subprocess,sys
import optparse
import re
import glob
import os
import shelve
#my function
from anfuncs import * 
import anfuncs
#Begin
if len(sys.argv)>1:# no argument-read from database lofl.s, 1 or more - read from real files
  allfiles = []
  #os.chdir('puregb') #
  fileList1 = glob.glob('puregb/*.out*')
  allfiles.extend(fileList1)
  fileList1 = glob.glob('C/*.out*')
  allfiles.extend(fileList1)
  fileList1 = glob.glob('N/*.out*')
  allfiles.extend(fileList1)
  fileList1 = glob.glob('O/*.out*')
  allfiles.extend(fileList1)
  fileList1 = glob.glob('Cbigz/*.out*')
  allfiles.extend(fileList1)
  fileList1 = glob.glob('Cbigy/*.out*')
  allfiles.extend(fileList1)
  fileList1 = glob.glob('noimp/*.out*')
  allfiles.extend(fileList1)
  fileList1 = glob.glob('Si/*.out*')
  allfiles.extend(fileList1)
  fileList1 = glob.glob('H/*.out*')
  allfiles.extend(fileList1)

  #Determing  number of files
  nfiles=len(allfiles)


  dig = re.compile('\d+') # make object w, mathed only with digits 
  j = 0; i = 0;n_st = 0;i_st = 0
  lofl = {};  lofse = [];  loft = []
  for dirandname in allfiles:
    nn = len(dirandname.split('/'))-1 # index of file name (always last) in list dirandname.split('/')
    name = dirandname.split('/')[nn] # name of file without dir 
    s1 = name.split(".") #split string by words using . as delimiter
    a1 = dig.search(s1[1]).group() # extract from s1[1] first group of digits
    if not s1[0].isdigit() or not a1.isdigit : #Check format of names
        print "Warning, file",  dirandname,"has bad name!, continue others"
        continue        
    it = int(s1[0]) #number of type
    ise = int(a1) #number of set
    print "it= ",it, " ise= ",ise
    if(not ise in lofse): 
        lofse.append(ise)
        i = i+1
        
    if(not it in loft): 
        loft.append(it)
        j = j+1
        
#    numit[it] = j
#    numise[ise] = j
    lofl[(it,ise)] = CalcResults()
    lofl[(it,ise)].iset=ise
    lofl[(it,ise)].itype=it
    lofl[(it,ise)].name = name
    lofl[(it,ise)].dirandname=dirandname
    lofl[(it,ise)].read_all(dirandname) 
    i_st = i_st + 1
#    lofl[(it,ise)].i_st(i_st)
  n_st = i_st
#Sort lists, Determing number of sets and types
  loft.sort()
  lofse.sort()
  d = shelve.open('lofl.s')
  d['lofl'] = lofl
  d['loft'] = loft
  d['lofse'] = lofse
  d.close()

d = shelve.open('lofl.s')
lofl = d['lofl']
loft = d['loft'] 
lofse = d['lofse']
d.close()
ntypes = len(loft)
nsets = len(lofse)
# num=it*nsets+ise
print "\nNumber of sets =",nsets, ",Some sets and types may be not useable, see class"
print "Number of types =",ntypes
#print "\n\nNew determing based on classes!!!"
#print "Energy charactericstics for every set between ideal Ti cell with 56 atoms and all other cells\n"



#force_analysis(lofl[(0,502)])
#print lofl[(601,502)].fcart
#force_analysis(lofl[(601,502)])
#force_analysis(lofl[(602,502)])
#force_analysis(lofl[(603,502)])
#write_md(lofl[(0,502)])

#sys.exit()



os.chdir('geoout')

#lofl[(606,502)].ietotal[len(lofl[(606,502)].ietotal)-1]=-3173.03395336953 #from 606.504

for ise in [502,500]: #lofse:
    print "SET=",ise
    for it in loft:
        try: 
            lofl[(it,ise)]
        except KeyError:
            continue
        if (lofl[(it,ise)].useable == 0): continue 
        cur = lofl[(it,ise)]
        print "for file %s: " % (cur.name) 
        # Calculate dE with E of Ti2 cell from "sravnenie faz s 2 atomami ... .xls"" file
        #Eti2 = -117.31019532 # Energy of two titanium atoms (primitive cell)
        #dE = (lofl[1][ise].etotal / lofl[1][ise].nznucl[0]) - (Eti2 * 0.5)
        
        #print "((Eti56-28*Eti2)/56) = dE = %g meV" % (dE*27.21138386*1000) #end of test-OK
        
        # Calculating energy of grain boundary including segregations
        # 1.Calculate energy of dope(impurity) in intersitials
        # Name conventions for types: 1) 600 for carbon in octa, 601-699 for other sites of carbon
        # 2) 700 for nitrogen in octa
        # 3) 800 for oxygen in octa  

        zcur = it/100
        if  cur.ntypat == 2:
            zcur = cur.znucl[1] #impurity charge of current structure set
            #if zcur not == lofl[zcur*100][ise].znucl[1]: # 600 for Carbon in octahedral site of ideal lattice  
            #    print "Warning! the name of file="lofl[zcur][ise].name," is not consistent with \
            #    it internal data, name of file must begin with znucl[1] = ",znucl[1]
            #    os._exit(1)

#Trying to find bulk cell of current ise and make Ebulk
        try:
            bk = lofl[(1,ise)]
            Ebulk = bk.etotal/bk.nznucl[0]                
        except KeyError:
            print 'No bulk energy for current set. Using from "sravnenie faz s 2 atomami ... .xls"'
#            bk = lofl[(1,502)]
            Ebulk = -58.65509766; #bk.etotal/bk.nznucl[0] 

#Trying to find pure gb cell  of current ise and make Ecellgb
        try:
            pgb = lofl[(0,ise)]
            Ecellgb = pgb.ietotal[len(pgb.ietotal)-1]
        except KeyError:
            print 'No energy of cell with pure GB. Setting to zero'
            Ecellgb = 0
          
#Trying to find bulk cell with impurity in octahedral site of current ise and make Eoct
                #Eoct = 5.8050671 # Ha, from  xls  for some time
        try:
            octbk = lofl[(zcur*100,ise)]
            Eoct = octbk.ietotal[len(octbk.ietotal)-1] - bk.etotal
        except (KeyError,NameError):
            print 'No oct energy for current set. Using from "sravnenie faz s 2 atomami ... .xls"'
            Eoct = -5.8050671

#Trying to find gb cell with impurity in 601 site of current ise and make  Eref
        try:
            irefgb = lofl[(zcur*100 + 1,ise)]
            Eref = irefgb.ietotal[len(irefgb.ietotal)-1] #Last full Energy of cell with atom almost between \
                           #grains (~11.8 Bhor) 
        except (KeyError,NameError):     
            print 'No energy for 601 ref state for current set. Setting to zero'       
            Eref = 0            

        if cur.ntypat == 1:
            Eoct = 0 #Because in current lofl there are no impurities
            zcur = 0 #Means no segregation exists

        
        #2.Calculate grain boundary energy with account of impurities
        #grain boundary is perpendicular to y axis
        print "z of impurity = ", zcur
        print "etotal every MD step! Eoct = ", Eoct
        w6 = len(cur.ietotal) - 3 #number of steps to show
        print cur.ietotal[0]
        
        for i in range(w6,len(cur.ietotal)):
            if i < 0 : continue
            cur.Egb = (cur.ietotal[i] - ( Ebulk * cur.nznucl[0]) - Eoct) / (2 * cur.sqrxz) * 97.17369323 *16 # J/M^2
            print "MD step = %i, Egb = %g J/M^2  Ediffmd = %g meV" % (i, cur.Egb,cur.Ediffmd[i]*27.21138386*1000)
#            if it == 0 : 
#                Egbpure = Egb
            if  cur.ntypat == 2 or 1:
                cur.Eseg = (cur.ietotal[i] - Eref) * 27.21138386 # eV
#                Eseg2 = (Egb - Egbpure) * cur.sqrxz * 2 #Calculating Eseg from Egb  #no sense, because impurity exist only on single gb
                Eseg3 = cur.ietotal[i] - Ecellgb  - Eoct #Must be the same as Eseg2 computed directly
                print "Eseg = %g eV, Eseg3 = %g eV" % (cur.Eseg,Eseg3*27.21138386 )

            #3Print etotal, max, rms force on atoms at the end of each MD step                  
            print cur.fcart[i],"\n"

#        if 1:#cur.calc_end: #Print only if calculation normaly completed
#            Egb = (cur.etotal - (Ebulk * cur.nznucl[0]) - Eoct) / (2 * cur.sqrxz) # 2 because \
        
#            print "End: znucl of segreg=%i, Egb=%g eV/A^2, Eoct=%g eV" % (zcur, Egb*97.17369323, Eoct*27.21138386)     
#            print cur.time,"\n"


#Make geo
for ise in lofse:
    for it in loft:
        try: 
            st = lofl[(it,ise)]
        except KeyError:
            continue
        if (st.useable == 0): continue 
        write_geo(st)

#out_for_paper(lofl,[106,102,103,104,105,101,606,602,603,604,605,601,706,702,703,704,705,701,806,802,803,804,805,801,1406,1402,1403,1404,1405,1401],[502,])

#out_for_paper(lofl,[0,606,602,603,604,601,706,702,703,704,701,806,802,803,804,801],[502,])


#addit_processing(lofl)








shift_analys(lofl[(601,702)])
#write_geo(lofl[(1401,502)])
#write_geo(lofl[(101,502)])


#write_geo(lofl[(1,504)])
#write_geo(lofl[(0,5022)])

#write_geo(lofl[(601,43)])
#write_geo(lofl[(602,43)])
#write_geo(lofl[(0,4)])    
#This coordinates made by geo_cell program!!!

#log
#19.07.2012 
#in class CalcResults command for xcart was changed. Now it use the last (xcart) "there are in .out for each MD step"
#!Pay attention that in construct.py the old version for xcart command is used, in which "xcart " is searched.
#!class CalcResults can not correctly read the last rprimd (rprim and acell) in the case
#of geometry relaxation and uncomplete calculation.
#xcart changed by xcart_str. added xcart[(i,j)]



