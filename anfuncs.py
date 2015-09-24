""" 
Include:
1. runBash(cmd)
2. CalcResults
3. interstitial()
4. out_for_paper()
5. shift_analys(st)
6. write_geo(st)
"""

import subprocess
import optparse
import re
import glob
import os
import math
import sys
import colorsys

#This function takes Bash commands and returns them
def runBash(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    out = p.stdout.read().strip()
    return out  #This is the stdout from the shell command


def write_md(st):
    f = open(st.name+".md",'w')   
    for i in 0, 1, 2:
        f.write('%f %f %f\n'%(st.rprimd.x[i],st.rprimd.y[i],st.rprimd.z[i]))
    f.write(str(st.natom)+'\n')
    f.write(st.typat.sdig+'\n')
    for i in range(st.natom):
#        if int(st.typat.v[i]) == 1:
        f.write(str(st.xcart[(i,0)]*0.529177)+" "+str(st.xcart[(i,1)]*0.529177)+" "+str(st.xcart[(i,2)]*0.529177)+"\n")

    for i in range(st.natom):
        f.write("47.88 ")
    f.write("\n")
    f.close()



def addit_processing(lofl):
    #print loft
    #out_for_paper(lofl,loft,[4])
    print "\n"
#shift_analys(lofl[(7051,4)])
    print "\n"
    shift_analys(lofl[(11,4)])
    print "\n"
    print lofl[(2,4)].xcart_str



    lofl[(0,532)].etotal = -6334.434391849

    muCdimond = -5.74110
    #For cell size convergence check 
    ise = 502
    it = 601
    itref = 0
    norm = lofl[(it,ise)].nznucl[0]*1.0/lofl[(itref,ise)].nznucl[0]
    Eoct = 0#(lofl[(600,ise)].etotal - lofl[(1,ise)].etotal)
    gbsqr = lofl[(it,ise)].acell[0] * lofl[(it,ise)].acell[2]
    Eseg1 = (lofl[(it,ise)].etotal - norm*lofl[(itref,ise)].etotal - Eoct)/gbsqr/2
    ise = 532
    norm = lofl[(it,ise)].nznucl[0]*1.0/lofl[(itref,ise)].nznucl[0]
    Eoct2 = (lofl[(600,ise)].etotal - lofl[(itref,ise)].etotal)
    gbsqr = lofl[(it,ise)].acell[0] * lofl[(it,ise)].acell[2]
    Eseg2 = (lofl[(it,ise)].etotal - norm*lofl[(itref,ise)].etotal - Eoct)/gbsqr/2


    print "Eoct2=",Eoct2
    print "Egb with carbon 601: 502, 532 =",Eseg1*97.17369323*16," ",Eseg2*97.17369323*16
    #print "Eseg 601: 502, 532 =",Eseg1*27.2," ",Eseg2*27.2
    #Egb without carbon 601: 502, 532 = 0.75293321353   0.968402913629
    print "fcart 601: 502, 532 =",lofl[(601,502)].fcart[len(lofl[(601,502)].fcart)-1], lofl[(601,532)].fcart[len(lofl[(601,532)].fcart)-1]


    #print "dif Eoct,eV zcell 502, 532 =",(Eoct-Eoct2)*27.2
    print lofl[(601,532)].acell[1]/lofl[(601,502)].acell[1]



#Convergence by y
    ise=54
    Eseg3 = lofl[(611,ise)].ietotal[len(lofl[(611,ise)].ietotal)-1] - lofl[(617,ise)].ietotal[len(lofl[(617,ise)].ietotal)-1]
    print "ise = ",ise
    print "diff 611-617=", Eseg3*27.2*1000
    print "fcart: 611, 617 =",lofl[(611,ise)].fcart[len(lofl[(611,ise)].fcart)-1], lofl[(617,ise)].fcart[len(lofl[(617,ise)].fcart)-1]
    
    Ex = lofl[(10,ise)].etotal/lofl[(10,ise)].natom - lofl[(0,502)].etotal/lofl[(0,502)].natom
    print "diff 10.ise-0.502=", Ex*27.2*1000

    it=10
    norm = lofl[(it,ise)].nznucl[0]*1.0/lofl[(1,502)].nznucl[0]
    gbsqr = lofl[(it,ise)].acell[0] * lofl[(it,ise)].acell[2]
    Egb = (lofl[(it,ise)].etotal - norm*lofl[(1,502)].etotal)/gbsqr/2
    print "Egb(10.ise)=", Egb*97.17369323*16,'J/m^2' 




#Analys of excess volume
    vTi2 = 234.082 #bohr^3 from tihcp1.5.out

#gb = lofl[(0,5022)]
    gb = lofl[(0,504)]

    bk = lofl[(1,504)]
#bk = lofl[(1,502)]
    gb.s = gb.acell[0]*gb.acell[2]
    gb.e = ( gb.v - ( gb.natom * (bk.v / bk.natom) ) )/gb.s/2*0.529177
    print "Excess volume of ",gb.name,' is ',gb.e , 'Angs '
#Estimation of influence of cell sizes
    print "Dif of volume Ti56-Ti2 on atom",(bk.v/bk.natom - vTi2/2)*gb.natom/gb.s/2*0.529177, ' Angs' # only 0.003 Angstrom, which means
#that 100 times smaller than typical excess volume of GB (0.3 A), but for 1.502 the difference is 0.03 Angstrom, which is quiet large,
#but still 10 times smaller; using vTi2 is good approch and there is no need to use 1.502






    #Show grain boundary regions A
    ise=502
    f = open('gb_structureA.xyz','w')
    f.write('40\n')
    f.write('gb_structure\n')
    xcart={}
    yshift = 0
    for it in 0,601,603,606:
        st=lofl[(it,ise)]

        for i in 5,40,49,24,14,4,13,23,48,39:
            i=i-1
            f.write('Ti '+str(st.xcart[(i,0)]*0.529177)+" "+str(st.xcart[(i,1)]*0.529177 + yshift)+" "+str(st.xcart[(i,2)]*0.529177)+'\n')
        yshift+=10
    f.close()
    #Show grain boundary regions B
    f = open('gb_structureB.xyz','w')
    f.write('36\n')
    f.write('gb_structure\n')
    xcart={}
    yshift = 0
    for it in 601,603,606:
        st=lofl[(it,ise)]

        for i in 8,16,26,25,15,6,7:
            i=i-1
            f.write('Ti '+str(st.xcart[(i,0)]*0.529177)+" "+str(st.xcart[(i,1)]*0.529177 + yshift)+" "+str(st.xcart[(i,2)]*0.529177)+'\n')
        for i in 51,50,41,34,33:
            i=i-1
            f.write('Ti '+str(st.xcart[(i,0)]*0.529177)+" "+str(st.xcart[(i,1)]*0.529177+yshift-st.acell[1]*0.529177)+" "+str(st.xcart[(i,2)]*0.529177)+'\n')

        yshift+=10
    f.close()

    #Show omega phase at grain boundary region B
    f = open('gb_structureBshowOmega.xyz','w')
    f.write('19\n')
    f.write('gb_structure\n')
    xcart={}
    yshift = 0
    st=lofl[(0,502)]
    copy(st,1)
    copy(st,2)
    copy(st,3)
    copy(st,3)
    for i in 816,600,483,482,566,781,698,772,565,556,340,771,555,548,762,546,466,689,473:
        i=i-1
        f.write('Ti '+str(st.xcart[(i,0)]*0.529177)+" "+str(st.xcart[(i,1)]*0.529177 + yshift)+" "+str(st.xcart[(i,2)]*0.529177)+'\n')
    f.close()

    #Show octahedral sites around carbon
    f = open('octahedral601_603.xyz','w')
    f.write('12\n')
    f.write('gb_structure\n')
    xcart={}
    yshift = 0

    for it in 601,:#701,801:
        st=lofl[(it,ise)]
        for i in 1,2,9:
            i-=1
            f.write('Ti '+str((st.xcart[(i,0)])*0.529177)+" "+str(st.xcart[(i,1)]*0.529177+yshift)+" "+str(st.xcart[(i,2)]*0.529177)+'\n')
        for i in 20,28,27:
            i-=1
            f.write('Ti '+str((st.xcart[(i,0)]-st.acell[0])*0.529177)+" "+str(st.xcart[(i,1)]*0.529177+yshift)+" "+str(st.xcart[(i,2)]*0.529177)+'\n')

    yshift+=10

    for it in 603,:#703,803:
        st=lofl[(it,ise)]
        for i in 21,23,29,22,24,14:
            i-=1
            f.write('Ti '+str(st.xcart[(i,0)]*0.529177)+" "+str(st.xcart[(i,1)]*0.529177 + yshift)+" "+str(st.xcart[(i,2)]*0.529177)+'\n')


    f.close()






def copy(st,axis):
    a = 0; b = 0; c = 0;
    if axis == 1: a = st.acell[0]
    if axis == 2: b = st.acell[1]
    if axis == 3: c = st.acell[2]
    for i in range(st.natom):
        st.xcart[(i+st.natom,0)] = st.xcart[(i,0)] + a
        st.xcart[(i+st.natom,1)] = st.xcart[(i,1)] + b 
        st.xcart[(i+st.natom,2)] = st.xcart[(i,2)] + c
        st.typat.sv.append(st.typat.sv[i])      
    st.acell[0] = st.acell[0] + a
    st.acell[1] = st.acell[1] + b    
    st.acell[2] = st.acell[2] + c 
    st.natom = 2*st.natom 

def force_analysis(st):
    max_r = max(st.fcart_end.r)
    h = [i/max_r/360.*300 for i in st.fcart_end.r]
#all values has ranges in 0-1
    s = 1; v = 1;
    f = open(st.name+".jmol",'w')
    f.write('load "visualisation/segreg_abinit/'+st.name+'.xyz"'+"\n")
    string = """
    set boundbox 2
    background white
    rotate z 90
    select all
  
    """
    f.write(string)
    r=[];g=[];b=[];
    for i in range(len(h)):
        rgb = colorsys.hsv_to_rgb(h[i], s, v)
        r.append(int(rgb[0] * 255))
        g.append(int(rgb[1] * 255))
        b.append(int(rgb[2] * 255))

        f.write('select Ti'+str(i+1)+"\n")
        f.write('color [%i,%i,%i]\n'%(r[i],g[i],b[i]))

    f.close()





comdic = {} #dictionary of commands




class DataStuff:
    """Parent Class For reading data and manipulate"""
    def do(self,command):
        self.s = runBash(command)
        ch = re.compile('[a-zA-Z]+') # make object , matched anything other than digit
        #find first entrance of any word
        try:
            w = ch.match(self.s).group()
            self.sdig = self.s.replace(w, '')#delete
            self.sv =  self.sdig.split() #values but like strings
            self.v = []
            for i in range(len(self.sv)):
                self.v.append(float(self.sv[i]))  
        except AttributeError:
            self.sdig = self.s
                
        
class CartArray(DataStuff):
    "Objects of this class will storage lists of 3 cartesian components, like list of atoms coordinates "
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
        self.r = []

    def do(self,command):
        DataStuff.do(self,command)
        if not self.sdig == '':
            lines = self.sdig.splitlines() #divide xcart string by lines
            for i in range(len(lines)):
                w = lines[i].split() #divide each line on components 
                if not len(w) == 3: 
                    print "Too many components, check sdig attribute - allowed only 3"
#                print self.s
                
                
                x = float(w[0]); y = float(w[1]); z = float(w[2]);
                self.x.append(x)
                self.y.append(y)
                self.z.append(z)
                self.r.append( math.sqrt(x*x + y*y + z*z) )

class CalcResults:

    def read_all(self,filename):
        self.useable = 0 
        self.calc_end = 0
        self.Eseg = 0
#read acell for two cases 
        command1="""grep -w "acell " """+filename+""" | tail -n 1"""
        command2="""grep -w -A 1 "(acell)" """+filename+""" | tail -n 1"""
        s1=runBash(command1)
        s2=runBash(command2)
#        print s1
#        print command2
#        print s2
        if s1=='':
            self.acell=[0,0,0]
            self.acell_str = 'no_acell'
            self.sqrxz = 0
            self.V = 0
        else:
            s = 'acell   ' + s2 
            if s2 == '': s = s1 

            s.strip()
            s_t = s.split()
            self.acell=[float(s_t[1]),float(s_t[2]),float(s_t[3])] 
            self.acell_str = s 
            self.sqrxz = self.acell[0] * self.acell[2] #S = X * Z, assumed that \
            self.v = self.acell[0] * self.acell[0]* self.acell[2]  #Only for cubic cells

#Find ucvol
        command="""grep -w "(ucvol)" """+filename+""" | tail -n 1"""
        s1=runBash(command)
        if s1=='':
            self.ucvol=0
            self.ucvol_str = 'no_ucvol'
        else:
            s1.strip()
            s_t = s1.split('=')
            self.ucvol = float(s_t[1])
            self.ucvol_str = s1



        command="""grep -w "ntypat " """+filename+""" | tail -n 1"""
        s1=runBash(command)
        self.ntypat_str = s1
        if s1 == '':
            self.ntypat = 0
        else:
            s3=s1.split()
            self.ntypat=int(s3[1]) 
            if(self.ntypat > 2):
                print "Program doesnt support more than two types of atoms!"
                os._exit(1)

        command="""grep -w "znucl " """+filename+""" | tail -n 1"""
#        s4 = ["", 0, 0]
        s1 = runBash(command)
        if s1 == '':
            print "no znucl, not useable"
            return
        s4 = s1.split() 
#        print s1
#        print self.ntypat
        self.znucl = []
        self.znucl_str = s1
        for i in range(self.ntypat):
            self.znucl.append(int(float(s4[i+1])))#i+1 because at i=0 lie word znucl
            #print s4
            #print self.znucl
          
        #Read total number of atoms and number of atoms of type znucl[0], which is assumed to be natom-1 if ntypat == 2
        command="""grep -w -m 1 "natom  " """+filename
        s1=runBash(command)
        self.natom_str = s1
        if s1=='':
            self.natom = 0
             
            print "!Caution In filename",filename,"not found natom! set to zero. It is very likly that other parameters was not \
            found too, Calculation completly not useable!!!"
        else:
            
            self.useable = 1
            self.nznucl = [] #Numbers of atoms by their types, for example 6 Ti atoms, 3 C atoms... 
            self.natom=int(s1.split()[1]) 

            if(self.ntypat == 1):
                self.nznucl.append(self.natom) 
            else:
                self.nznucl.append(self.natom-1)
                self.nznucl.append(1)

        command="""grep -w -m 1 -B 1 "+Overall" """+filename
        s1=runBash(command)
        if s1=='':
            self.time="Caution! Calculation not completed!!!"
            print "!In filename",filename,"Calculation not completed,!"
        else:
            self.time=s1 
            self.calc_end = 1
        
        #Find fcart, etotal and Absolute differense between etotals every MD   
        command="""grep -w "(fcart)" """+filename


        s1=runBash(command)
        if s1=='':
            print "Calculation has not fcart !!!"
            self.fcart = []
            self.fcart.append('No forces (fcart) !')

        else:
            self.fcart = s1.splitlines() 
  
#find fcart_end        
        nMDOK = 1
        if len(self.fcart) > 1:
            nMDOK = len(self.fcart) - 1
        comdic['fcart_end'] = 'grep -w -A ' + str(self.natom) + ' -m ' + str(nMDOK) + \
            ' "(fcart)" ' + filename + ' | tail -n ' + str(self.natom) + ';'
        print comdic['fcart_end']
        self.fcart_end = CartArray(); 
        self.fcart_end.do(comdic['fcart_end'])

#        sys.exit()  
  
  
  
             

        command="""grep -w "etotal " """+filename+""" | (read a b; echo $b)"""
        s1=runBash(command)
        if s1=='':
            self.etotal=0
            print "!In filename",filename,"not found etotal! set to zero"
        else:
            self.etotal=float(s1) 

# find etotal every md step
        command="""grep -w "(etotal)" """+filename
        s1=runBash(command)
        if s1=='':
            print "Calculation has no (etotal) !!!"
            self.ietotal = []; self.ietotal.append(self.etotal)
        else:
            s5 = s1.splitlines()
            self.ietotal = []
            for i in range(len(s5)): self.ietotal.append(float(s5[i].split('=')[1]))
 
 
        
        command="""grep -w "Absolute" """+filename
        s1=runBash(command)
        if s1=='':
            print "Calculation has not (etotal) !!!"
            self.Ediffmd = [0]
        else:
            s6 = s1.splitlines()
            self.Ediffmd = []
            self.Ediffmd.append(self.ietotal[0])
            for i in range(len(s6)): self.Ediffmd.append(float(s6[i].split('=')[1]))

        #read xcart - if natom > 0        
        if self.natom > 0:
            command1 = 'grep -w -A ' + str(self.natom - 1) + \
            ' "xcart " ' + filename + ' | tail -n ' + str(self.natom) + '; '
            command2 = 'grep -w -A ' + str(self.natom) + \
            ' "(xcart)" ' + filename + ' | tail -n ' + str(self.natom) + ';'
            
#            print command1
#            print command2
            s1=runBash(command1)
            s2=runBash(command2)
            
            self.xcart_str = "xcart     " + s2
            if s2 == '':self.xcart_str = s1
            self.xcart = {}
            xcart_lines = []
            xcart_t = self.xcart_str.replace('xcart', '') #replace all xcart by nothing in string
            xcart_lines = xcart_t.splitlines() #divide xcart string by lines
            for i in range(self.natom):
                for j in 0, 1, 2:
                    self.xcart[(i,j)] = float(xcart_lines[i].split()[j]) #divide each line on components and transform to float
            
            
            if s2 == '':
                self.xcart_str = s1


#            
 #               print self.xcart
        else:
            self.xcart_str = "no xcart"         
            self.xcart = 0



#read typat
        command="""grep -w -A 2 "typat " """+filename+""" | tail -n 3"""
        if self.natom > 100:command = """grep -w -A 5 "typat " """+filename+""" | tail -n 6"""
        self.typat = DataStuff()
        self.typat.do(command)
        if self.typat.s=='':
            print "Calculation has no typat !!!"





        comdic['rprim'] = """grep -w -A 2 "rprim " """+filename+""" | tail -n 3"""
        self.rprim = CartArray()
        self.rprim.do(comdic['rprim'])
        self.rprim_str = self.rprim.s

        if self.rprim.sdig == '':
            self.rprim.x.append(1);self.rprim.y.append(0);self.rprim.z.append(0)
            self.rprim.x.append(0);self.rprim.y.append(1);self.rprim.z.append(0)
            self.rprim.x.append(0);self.rprim.y.append(0);self.rprim.z.append(1)

        self.rprimd = CartArray()
        for i in 0, 1, 2:
            self.rprimd.x.append(self.acell[i]*self.rprim.x[i]);
            self.rprimd.y.append(self.acell[i]*self.rprim.y[i]);
            self.rprimd.z.append(self.acell[i]*self.rprim.z[i])





def write_geo(st):
    f = open(st.name+".geo",'w')
    f.write(st.acell_str+"\n"+st.rprim.s+"\n"+st.xcart_str+"\n"+st.ntypat_str\
    +"\n"+st.typat.s+"\n"+st.natom_str\
    +"\n"+st.znucl_str)
    f.close()



def write_xyz(st):
    f = open(st.name+".xyz",'w')
    f.write(str(st.natom)+"\n")
    f.write(st.name+"\n")
    for i in range(st.natom):
        if st.typat[i] == 1:
            f.write("Ti "+str(st.xcart[(i,0)]*0.529177)+" "+str(st.xcart[(i,1)]*0.529177)+" "+str(st.xcart[(i,2)]*0.529177)+"\n")

            
            
        if st.typat[i] == 2:
            f.write("C "+str(st.xcart[(i,0)]*0.529177)+" "+str(st.xcart[(i,1)]*0.529177)+" "+str(st.xcart[(i,2)]*0.529177)+"\n")  
    f.close()



def interstitial():
    #This coordinates made by geo_cell program!!!
    #Making geometry files for cells with intersetials in different positions
    #14.2159 11.7811 2.18529 Atom Number 0   601
    #1.60525 22.7022 2.18529 Atom Number 1   602
    #11.0632 20.8821 2.18529 Atom Number 2   603
    #5.8088 19.0619 2.18529 Atom Number 3    604
    #15.2668 17.2417 2.18529 Atom Number 4   605
    #5.88496 25.4826 0 -coordinates of site #1 in the plane of GB  606 


    #This coordinate was used for octahedral site natom in cell without gb
    #5.8088 19.0619 2.18529 Atom Number 3   600

    #Coordinates of intersetials
    inters_cart = [0 for i in range(10)]
    inters_cart[1] = "\n14.2159 11.7811 2.18529" #Atom Number 0   601   11.8 Bhor from grain
    inters_cart[2] = "\n1.60525 22.7022 2.18529" #Atom Number 1   602   2.8 Bhor from grain
    inters_cart[3] = "\n11.0632 20.8821 2.18529" #Atom Number 2   603   4.6 Bhor from grain
    inters_cart[4] = "\n5.8088 19.0619 2.18529" #Atom Number 3    604   6.5 Bhor from grain
    inters_cart[5] = "\n15.2668 17.2417 2.18529" #Atom Number 4   605   8.3 Bhor from grain
    inters_cart[6] = "\n5.88496 25.4826 0" #-coordinates of site #1 in the plane of GB  606 


    #Making geo files for different positions
    st = lofl[(0,4)] #Maximum relaxed grain boundary without any impurities, Do not copying!!!!!
    #from 0.4yesgb.out 12 md step 
    xcart_base = """
xcart
-2.18498971344626E-01  8.18150142875686E+00  4.30092111809256E+00
  7.75342790634931E-01  1.38456280020084E+01  4.19513572181355E+00
  1.81711204493855E+00  1.95503965452007E+01  4.06316563504539E+00
 -5.50047342881612E-01  2.54359541945873E+01 -5.02744517094890E-02
  2.99734964189795E+00  2.54364440818155E+01  4.39766599808313E+00
  1.70092502582476E+00  2.96867846469700E-02 -3.02691271993559E-02
  1.53972502321383E+00  5.32313886659516E+00 -3.97586415135043E-02
  4.15627139259447E+00  3.02103001205849E+00  4.31797889280254E+00
  1.80472311272683E+00  1.08549074532435E+01 -1.21758507505547E-01
  4.87381702949080E+00  9.77301369743774E+00  4.37451261661295E+00
  2.74245947571436E+00  1.65074020373177E+01 -2.94481469405892E-01
  5.96185343618244E+00  1.56582764828773E+01  4.15045215745461E+00
  3.73676288888429E+00  2.22085556860016E+01 -2.67294256632624E-01
  7.01588784687555E+00  2.12440716217597E+01  3.67080576489632E+00
  6.63262441889292E+00  2.66612231316594E+00 -8.38070307324486E-02
  8.99628963077618E+00  2.75842590931577E-02  4.36166963345155E+00
  6.16624371462286E+00  7.67229058621607E+00  5.73732869392318E-02
  9.07411060879979E+00  5.62899470218345E+00  4.38143165910889E+00
  6.95675936064015E+00  1.27742859497289E+01 -1.81795577323557E-01
  1.01332709168137E+01  1.14984678442532E+01  4.25750721572554E+00
  7.97974608988828E+00  1.83047935684179E+01 -4.41905731982096E-01
  1.11891391294477E+01  1.72576691535552E+01  4.11621383892166E+00
  8.79809515370853E+00  2.40158731700907E+01 -1.62275492485192E+00
  1.22824114440917E+01  2.27828097453895E+01  3.91705468208765E+00
  1.14759563410213E+01  2.75520703783484E+00 -1.16839777315812E-02
  1.38934312079851E+01  2.72241138275434E+00  4.36824874355055E+00
  1.13423425098078E+01  8.40974480128966E+00 -7.89984902500356E-02
  1.21436755439857E+01  1.42837777317396E+01 -2.31463775083823E-01
  1.31976473644886E+01  1.98252264394924E+01 -3.33164233543386E-01
 -5.45384640535953E-02  4.27368234764728E+01  4.42788752727697E+00
  9.39078579710003E-01  3.70467224478855E+01  4.53210488023706E+00
  1.97361873385470E+00  3.13535267865895E+01  4.70996254929467E+00
  1.59533700065463E+00  4.56192861787809E+01  1.47970754833347E-02
  4.19255976819220E+00  4.78707752788851E+01  4.39193902851551E+00
  1.93017214922937E+00  4.00625196160757E+01  8.80970368389759E-02
  5.01212764030189E+00  4.10510764377906E+01  4.45108257379112E+00
  2.96427484779129E+00  3.44112387641458E+01  2.31538109923039E-01
  6.13216175751316E+00  3.52763258523192E+01  4.62379147833618E+00
  3.83193117663702E+00  2.85744248035306E+01  2.50222618987879E-01
  7.12156578675259E+00  2.96708547879906E+01  5.11430721050934E+00
  6.65380756745582E+00  4.82944528681627E+01 -7.37341228528807E-04
  6.23507358380167E+00  4.32894933848349E+01  2.17685316352953E-02
  9.15143295117319E+00  4.53005423269813E+01  4.38602813813631E+00
  7.12782653615017E+00  3.81648324083008E+01  1.66399415559369E-01
  1.03050019191890E+01  3.94010399511393E+01  4.51668607886961E+00
  8.13140394937588E+00  3.25680835087574E+01  4.34215313666389E-01
  1.13564369278969E+01  3.36375646543716E+01  4.62997429586444E+00
  8.89312718881160E+00  2.68989895145841E+01  1.61974013647065E+00
  1.23995080995420E+01  2.81827762818521E+01  4.84793239555703E+00
  1.15024975737144E+01  4.81728264006837E+01 -3.07239400870196E-02
  1.39338038054685E+01  4.82019353503103E+01  4.35897540449638E+00
  1.14406190202066E+01  4.25069625522842E+01  3.33223660712672E-02
  1.23141597718047E+01  3.66227426558026E+01  1.73923761193133E-01
  1.33715067301235E+01  3.11011482902399E+01  2.88463735264039E-01
"""

    st.ntypat_str = "ntypat 2"
    st.typat = st.typat + " 2"
    st.natom_str = "natom " + str(st.natom + 1) 


    #Additional code making out.geo not_relaxed files and for help in config file
    iset = 401
    for j in 6, 7, 8:
        st.znucl_str = "znucl 22 "+str(j)
        for i in 1, 2, 3, 4, 5, 6:
            st.name = str( (j * 100) + i ) + "." + str(iset) +".yesgb.out"
            st.xcart_str = xcart_base + inters_cart[i]
            #write_geo(st)
            runname = str( (j * 100) + i ) + "." + str(iset) +".yesgb"
            geoname = str( (j * 100) + i ) + "." + str(4) +".yesgb.out.geo"
    #        print "name=%s" % (runname)
    #        print "geo=%s" % (geoname)
    #        print "make_run"

    
    #tests
    #for i in 1, 2, 3, 4, 5, 6:
    #print float(inters_cart[i].split()[0])*0.529177," ", float(inters_cart[i].split()[1])*0.529177," ", float(inters_cart[i].split()[2])*0.529177














def shift_analys(st):
    init_gb = CalcResults()
    init_gb.natom = 54
    init_gb.xcart_str = """1.05087e-11 8.49421 4.37058
1.05089 13.9548 4.37058
2.10177 19.4153 4.37058
0 25.4826 0
3.15266 25.4826 4.37058
14.7124 50.9653 0
1.05089 5.46056 0
4.20354 4.85383 4.37058
2.10177 10.9211 0
5.25443 10.3144 4.37058
3.15266 16.3817 0
6.30531 15.775 4.37058
4.20354 21.8423 0
7.3562 21.2355 4.37058
5.25443 1.82019 0
8.40709 50.9653 4.37058
6.30531 7.28075 0
9.45797 6.67402 4.37058
7.3562 12.7413 0
10.5089 12.1346 4.37058
8.40709 18.2019 0
11.5597 17.5951 4.37058
9.45797 23.6624 0
12.6106 23.0557 4.37058
10.5089 3.64038 0
13.6615 3.03365 4.37058
11.5597 9.10094 0
12.6106 14.5615 0
13.6615 20.0221 0
1.05087e-11 42.471 4.37058
1.05089 37.0105 4.37058
2.10177 31.5499 4.37058
1.05089 45.5047 0
4.20354 46.1114 4.37058
2.10177 40.0441 0
5.25443 40.6509 4.37058
3.15266 34.5836 0
6.30531 35.1903 4.37058
4.20354 29.123 0
7.3562 29.7297 4.37058
5.25443 49.1451 0
6.30531 43.6845 0
9.45797 44.2912 4.37058
7.3562 38.2239 0
10.5089 38.8307 4.37058
8.40709 32.7634 0
11.5597 33.3701 4.37058
9.45797 27.3028 0
12.6106 27.9095 4.37058
10.5089 47.3249 0
13.6615 47.9316 4.37058
11.5597 41.8643 0
12.6106 36.4038 0
13.6615 30.9432 0
"""
    init_gb.acell = [14.7124, 50.9653, 8.74117];

#init_gb.xcart_db = {}
    xcart_lines = []
    init_gb.xcart = {}
#xcart_t = init_gb.xcart_str.replace('xcart', '') #replace all xcart by nothing in string
    xcart_lines = init_gb.xcart_str.splitlines() #divide xcart string by lines
    for i in range(init_gb.natom):
        for j in 0, 1, 2:
            init_gb.xcart[(i,j)] = float(xcart_lines[i].split()[j])
 
 
 

#Here begins calculation of shifts in x and z directions of one grain against other by atoms in centers of grains
    reg1 = {}
    reg2 = {}
    reghs = 6 #Half of size of regin of shift analysis
    k=0;l=0;
    for i in range(init_gb.natom):
        if init_gb.xcart[(i,1)]>init_gb.acell[1]/4.-reghs and init_gb.xcart[(i,1)]<init_gb.acell[1]/4.+reghs:
            reg1[(k,0)] = init_gb.xcart[(i,0)]
            reg1[(k,2)] = init_gb.xcart[(i,2)]
            reg1[(k,3)] = i #number of atom
            k=k+1
        if init_gb.xcart[(i,1)]>(3*init_gb.acell[1]/4.)-reghs and init_gb.xcart[(i,1)]<(3*init_gb.acell[1]/4.)+reghs:
            reg2[(l,0)] = init_gb.xcart[(i,0)]
            reg2[(l,2)] = init_gb.xcart[(i,2)]
            reg2[(l,3)] = i
            l=l+1
    nreg1 = k
    nreg2 = l
    dx1 =0;dz1=0;dx2=0;dz2=0;
    for i in range(nreg1):
        dx1 = dx1 + reg1[(i,0)] - st.xcart[(reg1[(i,3)],0)]
        dz1 = dz1 + reg1[(i,2)] - st.xcart[(reg1[(i,3)],2)]
    dx1 = dx1/nreg1
    dz1 = dz1/nreg1

    for i in range(nreg2):
        dx2 = dx2 + reg2[(i,0)] - st.xcart[(reg2[(i,3)],0)]
        dz2 = dz2 + reg2[(i,2)] - st.xcart[(reg2[(i,3)],2)]
    dx2 = dx2/nreg2
    dz2 = dz2/nreg2
     
    dx=dx1-dx2
    dz=dz1-dz2
#    print "nreg1 = ",nreg1," nreg2 = ",nreg2
#    print "dx1 = ",dx1," dz1 = ",dz1
#    print "dx2 = ",dx2," dz2 = ",dz2
#    print "dx = ",dx," dz = ",dz
#    print dx*0.529177,"&",dz*0.529177,
    print '%2.2f & %2.2f'%(dx*0.529177,dz*0.529177),
    return

#Calculate distance from interstitial to GB
def distance_to_gb(st):
    dist = 0.
    try:
        dist = st.acell[1] / 2 - st.xcart[(54,1)]
    except KeyError:
        print "No st.xcart[(54,1)]"
#        print "distance from GB = ",dist*0.529177," A"
#    print dist*0.529177
    print '%2.2f'%(dist*0.529177),
    
    
    
    


def out_for_paper(lofl,loft,lofse):
    #Print dx and dz
    print "Filename, dx, dz, distance to Gb, Eseg, Eseg_mechanical, Eseg_chemical, dV; All in Angstroms and eV"
    for ise in lofse:
        for it in loft:
            try:
                cur = lofl[(it,ise)]
            except KeyError:
                continue
            try:
                Esegmech = lofl[(it,500)].Eseg
            except KeyError:
                Esegmech = 0

            lofl[(0,502)].znucl.append(0)
            if (cur.useable == 0): continue 
#            write_geo(cur)
            if cur.natom == 55 or 54:
                print cur.name,"&",;shift_analys(cur);print "&",;distance_to_gb(cur)
                print '& %2.2f & %2.2f & %2.2f & %2.2f'%(cur.Eseg, Esegmech, cur.Eseg - Esegmech, \
                (cur.ucvol - lofl[(cur.znucl[1]*100 + 1, 502)].ucvol) * 0.148184534 ) # in Angstroms
#    ise=502
#    for it in 606, 602, 603, 604, 605, 601:
#        distance_to_gb(cur); 
#        print ' %2.3f '%(cur.Eseg),;
#        print ' %2.3f '%(lofl[(it+100,ise)].Eseg),;
#        print ' %2.3f '%(lofl[(it+200,ise)].Eseg);
        
        
        
