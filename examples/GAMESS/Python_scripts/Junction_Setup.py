import subprocess, string, numpy

#This script is for setting up basic RUQT calculations from Gamess data. It is outdated so use with caution. JunctionMod and/or the Maple scripts are better choices.

print("Junction transport questionnaire!")

print("\n\nWhat is the name of your file? No file extensions:\nEXAMPLE: Benzene_jxn")
name = input("")

print("\nDo you have files named "+name+"T2, "+name+".gamess1, and "+name+".gamess2? [y/n]")
fileExists = input("")
if "y" in fileExists :
 print("Great, lets continue\n")
if "n" in fileExists :
 print("These files are needed in this directory for the script to work\n")
 quit()

print("\nWhat type of calculation will you perform?\nEXAMPLES: current, transmission")
while True :
 Calc_Type = input("")
 if Calc_Type in ["current", "transmission"] :
  break
 else :
  print("\nWhat you have entered was not valid. Please try again")

print("\nWhat software would you like to use?\nEXAMPLES: gamess, qchem, pyscf")
while True :
 qchem = input("")
 if qchem in ["gamess","qchem","pyscf"] :
  break
 else :
  print("\nWhat you have entered was not valid. Please try again")

print("\nWould you like double excitations?\nEXAMPLES: T, F")
while True :
 doubles = input("")
 if doubles in ["T","F"] :
  break
 else :
  print("\nWhat you have entered was not valid. Please try again")

print("\nWhat functional do you prefer?\nEXAMPLES: rdm, hf, cisd")
while True :
 functional = input("")
 if functional in ["rdm","hf","cisd"] :
  break
 else :
  print("\nWhat you have entered was not valid. Please try again")

print("\nHow many atoms for left electrode?")
elenum_l = int(input(""))
print("\nHow many atoms for right electrode?")
elenum_r = int(input(""))
#print("\nHow many atoms in the device region?")
#devnum = int(input(""))
devnum = 0 # this can be determined, makes life easier :)
print("\nThere are several more settings that are typically the same.\n Would you like to assume default values? y/n\nNote: 'y' is currently the only option")
isDefault = input("")

Electrode_Type = "Metal_WBL"
Fermi_enl = -5.5
Fermi_enr = -5.5
localden_fermi_l = 0.07
localden_fermi_r = 0.07
energy_start = -15
energy_end = 5
delta_en = 0.1
volt_start = 0
volt_end = 1
delta_volt = 0.1
KT = 0.00000000107

if isDefault == "n" :
 print("Sorry, but there is currently no option to setting non-default values here.\n You will have to rewrite it in the main jxn file")

print("Beginning search for remaining values:")
norb = 0
natomic = 0
numact = 0
numocc = 0
numvirt = 0
corr_ener = 0
with open(name + ".gamess2") as filesearch: #Edit for any file
 for line in filesearch:
  if "NUMBER OF CARTESIAN GAUSSIAN" in line:
   norb = int(line.split()[-1])
   natomic = norb 
   numact = norb 
  if "NUMBER OF OCCUPIED ORBITALS (ALPHA) KEPT" in line:
   numocc = int(line.split()[-1])
   numvirt = norb - numocc
   break
if functional == "rdm":
 with open(name + ".rdm") as filesearch:
  for line in filesearch:
   if "P2RDM CORRELATION ENERGY:" in line:
    corr_ener = float(line.split()[-1])
   break
#Determination of end_iter
line = ""
new_iter = ""
old_iter = ""
negf_mo = open(name + '.gamess1','r')
while 'DENSITY CONVERGED' not in line:
 line = negf_mo.readline()
 if 'END OF DEBUG' in line:
  old_iter = new_iter 
  new_iter = negf_mo.readline()
 if 'FINAL RHF ENERGY' in line:
  break
negf_mo.close()
end_iter = old_iter.split()[0] + " " + old_iter.split()[1]
#
size_c = 0
size_l = 0
size_r = 0
#Finding orbitals in device region and electrodes' regions ^ 3 variables above
with open(name + ".gamess2") as filesearch:
  numatoms = 0
  skipper = 0 #for skipping text after POP
  orbcount = 0
  curratom = 0
  c_isFound = False
  l_isFound = False
  die = False
  for line in filesearch:
    if 'TOTAL NUMBER OF ATOMS' in line:
      numatoms = int(line.split()[5])
      devnum = numatoms - elenum_l - elenum_r #devnum found here!!!
      continue
    if 'POPULATIONS IN EACH AO' in line:
      skipper = 1
      continue
    if skipper == 1 and line.split()[0] in 'MULLIKEN':
      skipper = 2
      continue
    linelist = line.split()
    if skipper == 2:
        ###left off here
      if len(linelist) == 6:
        curratom = int(linelist[2])
      if len(linelist) == 5:
        temp = linelist[1]
        curratom = int(temp[-2:])
      orbcount = orbcount + 1
      if curratom == (elenum_l+1) and l_isFound == False:
        size_l = orbcount - 1
        orbcount = 0
        l_isFound = True
      if curratom == (devnum+elenum_l+1) and c_isFound == False:
        size_c = orbcount
        orbcount = 0
        c_isFound = True
        size_r = norb - size_c - size_l
        die = True
    if die == True:
      break

#

negf_inp = open(name,'w')
negf_inp.write(Calc_Type +  "\n")
negf_inp.write(Electrode_Type +  "\n")
negf_inp.write("{0}".format(Fermi_enl) + "\n")
negf_inp.write("{0}".format(Fermi_enr) + "\n")
negf_inp.write("{0}".format(localden_fermi_l) + "\n")
negf_inp.write("{0}".format(localden_fermi_r) + "\n")
negf_inp.write("{0}".format(norb) + "\n")
negf_inp.write("{0}".format(natomic) + "\n")
negf_inp.write("{0}".format(numact) + "\n")
negf_inp.write("{0}".format(numocc) + "\n")
negf_inp.write("{0}".format(numvirt) + "\n")
negf_inp.write("{0}".format(size_c) + "\n")
negf_inp.write("{0}".format(size_l) + "\n")
negf_inp.write("{0}".format(size_r) + "\n")
negf_inp.write("{0}".format(energy_start) + "\n")
negf_inp.write("{0}".format(energy_end) + "\n")
negf_inp.write("{0}".format(delta_en) + "\n")
negf_inp.write("{0}".format(volt_start) + "\n")
negf_inp.write("{0}".format(volt_end) + "\n")
negf_inp.write("{0}".format(delta_volt) + "\n")
negf_inp.write("{0}".format(KT) + "\n")
negf_inp.write("{0}".format(qchem) + "\n")
negf_inp.write("{0}".format(doubles) + "\n")
negf_inp.write("{0}".format(functional) + "\n")
negf_inp.close()

negf_mo = open(name + '.gamess1','r')
#negf_htwo = open(name + '.Htwo','w')
#negf_smat = open(name + '.Smat','w')

negf_mo.seek(0)
index_col = int(norb/5)
index_col2 = numpy.float64(norb)/5-numpy.float64(index_col)
int_reman = numpy.int(index_col2)

index = numpy.zeros(10, dtype=numpy.int8, order='F')
temp  = numpy.zeros(10, dtype=numpy.float64, order='F')
temp2 = numpy.zeros(10, dtype=numpy.float64, order='F')
mo_ener = numpy.zeros(norb, dtype=numpy.float64, order = 'F')
mo_coeff = numpy.zeros((natomic,norb), dtype=numpy.float64, order='F')
#print len(lines)
#First we parse out the mo energies

line = negf_mo.readline()
while 'EIGENVECTORS' not in line:
 line = negf_mo.readline()

for y in range(0,3):
 line = negf_mo.readline()

index = line.split()
line2 = negf_mo.readline()
temp = line2.split()
#print index
#print temp
h = len(index)
for x in range(0,h):
 z = int(index[x])-1
 mo_ener[z] = temp[x]

if natomic > 5:
 for y in range(0,index_col):
  for x in range(0,natomic+3):
   line = negf_mo.readline()
  #for x in range(0,3):

  index = line.split()
  #print index
  line2 = negf_mo.readline()
  temp = line2.split()
  #print temp
  h = len(index)
  for x in range(0,h):
   z = int(index[x])-1
   mo_ener[z] = temp[x]

if index_col2 != 0:
 for y in range(0,1):
  line = negf_mo.readline()

 index = line.split()
 line2 = negf_mo.readline()
 temp = line.split()
 for x in range(0,int_reman):
  z = index[x]
  mo_ener[z] = temp[x]

#Now we get the mo coeff

#print mo_ener
negf_mo.seek(0)
while 'EIGENVECTORS' not in line:
 line = negf_mo.readline()

line = negf_mo.readline()
for y in range(0,index_col):
 line = negf_mo.readline()
 line = negf_mo.readline()
 if norb >= 4:
  index = line.split()
  #print index
  for w in range(0,2):
   line = negf_mo.readline()
  for x in range(0,natomic):
   line3 = negf_mo.readline()
   temp2 = line3.split()
   #print 't',x,temp2
   if 'X'==temp2[2] or 'Y'==temp2[2] or 'Z'==temp2[2] or 'S'==temp2[2]:
    w=3
   elif 'XX'==temp2[2] or 'YY'==temp2[2] or 'XY'==temp2[2] or 'XZ'==temp2[2]:
    w=3
   elif 'ZZ'==temp2[2] or 'YZ'==temp2[2]:
    w=3
   else:
    w=4
   h=len(index)
   for z in range(0,h):
    #print z,index[z],x,w,temp2[w]
    j = int(index[z]) - 1
    #print x,j#,w,temp2[2]
    mo_coeff[x,j] = temp2[w]
    w = w + 1

if index_col2 != 0:
 line2 = negf_mo.readline()
 line2 = negf_mo.readline()
 index = line2.split()
 h = len(index)
 for w in range(0,2):
  line = negf_mo.readline()
 for x in range(0,natomic):
  line3 = negf_mo.readline()
  temp2 = line3.split()
  #print 't',x,temp[2],temp[3]
  if 'X'==temp2[2] or 'Y'==temp2[2] or 'Z'==temp2[2] or 'S'==temp2[2]:
   w=3
  elif 'XX'==temp2[2] or 'YY'==temp2[2] or 'XY'==temp2[2] or 'XZ'==temp2[2]:
   w=3
  elif 'ZZ'==temp2[2] or 'YZ'==temp2[2]:
   w=3
  else:
   w=4
  for z in range(0,h):
   j = int(index[z]) - 1
   mo_coeff[x,j] = temp2[w]
   w = w + 1

#Print to file


negf_modat = open(name + '.mo_dat','w')
for x in range(0,natomic):
 for y in range(0,norb):
  #negf_modat.write("{0}".format(x+1)+" "+"{0}".format(y+1)+" "+"{0}".format(mo_coeff[x,y]) + "\n")
  negf_modat.write("{0}".format(mo_coeff[x,y]) + "\n")

for x in range(0,norb):
  #negf_htwo.write("{0},".format(x) + "{0},".format(y) +  "{0}".format(Htwo[x,y]) + "\n")
  negf_modat.write("{0}".format(mo_ener[x]) + "\n")

negf_modat.write("{0}".format(corr_ener) + "\n")

negf_modat.close()

## THIS IS THE SECOND FILE
negf_gamess1 = open(name + '.gamess1','r')

negf_gamess1.seek(0)
index_col = int(natomic/5)
index_col2 = numpy.float64(natomic)/5-numpy.float64(index_col)

index = numpy.zeros(10, dtype=numpy.int8, order='F')
temp  = numpy.zeros(10, dtype=numpy.float64, order='F')
Htwo = numpy.zeros((natomic,natomic), dtype=numpy.float64, order='F')
#print len(lines)

line = negf_gamess1.readline()
while end_iter not in line:
 line = negf_gamess1.readline()
#for line in negf_gamess1:
while 'TOTAL FOCK OPERATOR' not in line:
 line = negf_gamess1.readline()

for y in range(0,index_col):
 line = negf_gamess1.readline()
 line2 = negf_gamess1.readline()
 line = negf_gamess1.readline()
 if norb >= 4:
  index  =  line2.split( )
  #print index
  si =numpy.int16(index[0])-1
  #print si,natomic
  for x in range(si,natomic):
   line3 = negf_gamess1.readline()
   temp = line3.split()
   #print 't',x#,temp
   if 'X'==temp[2] or 'Y'==temp[2] or 'Z'==temp[2] or 'S'==temp[2]:
    w=3
   elif 'XX'==temp[2] or 'YY'==temp[2] or 'XY'==temp[2] or 'XZ'==temp[2]:
    w=3
   elif 'ZZ'==temp[2] or 'YZ'==temp[2]:
    w=3
   else:
    w=4
   n=len(temp)-w
   for z in range(0,n):
   #print z,index[z],x,w,temp
    j = int(index[z]) - 1
    k = int(temp[0])
    #print x+1,j+1,temp[w]
    Htwo[x,j] = temp[w]
    #if x != j:
    Htwo[j,x] = temp[w]
    w = w + 1

if index_col2 != 0:
 line = negf_gamess1.readline()
 line2 = negf_gamess1.readline()
 line = negf_gamess1.readline()
 index = line2.split()
 #h = len(index)
 #print index
 si =numpy.int16(index[0])-1
 for x in range(si,natomic):
  line3 = negf_gamess1.readline()
  temp = line3.split()
  #print 't',x,temp
  if 'X'==temp[2] or 'Y'==temp[2] or 'Z'==temp[2] or 'S'==temp[2]:
   w=3
  elif 'XX'==temp[2] or 'YY'==temp[2] or 'XY'==temp[2] or 'XZ'==temp[2]:
   w=3
  elif 'ZZ'==temp[2] or 'YZ'==temp[2]:
   w=3
  else:
   w=4
  n=len(temp)-w
  for z in range(0,n):
   j = int(index[z]) - 1
   k = int(temp[0])
   Htwo[x,j] = temp[w]
   #if x != j:
   Htwo[j,x] = temp[w]
   w = w + 1

negf_Htwo = open(name + '_Htwo','w')
for x in range(0,natomic):
 for y in range(0,natomic):
  #negf_Htwo.write("{0},".format(x+1) + "{0},".format(y+1) +  "{0}".format(Htwo[x,y]) + "\n")
  negf_Htwo.write("{0}".format(Htwo[x,y]) + "\n")

negf_Htwo.close()
negf_gamess1.close()

negf_gamess2 = open(name + '.gamess2','r')

index_col = int(natomic/5)
index_col2 = numpy.float64(natomic)/5-numpy.float64(index_col)

Smat = numpy.zeros((natomic,natomic), dtype=numpy.float64, order = 'F')

line = negf_gamess2.readline()
while 'OVERLAP MATRIX' not in line:
 line = negf_gamess2.readline()
#for line in negf_gamess1:
# if 'Overlap Matrix' in line:
#  break

for y in range(0,index_col):
 line = negf_gamess2.readline()
 line2 = negf_gamess2.readline()
 line = negf_gamess2.readline()
 if norb >= 4:
  index  =  line2.split( )
  #print index
  si =numpy.int16(index[0])-1
  for x in range(si,natomic):
   line3 = negf_gamess2.readline()
   temp = line3.split()
   #print 't',x,temp
   if 'X'==temp[2] or 'Y'==temp[2] or 'Z'==temp[2] or 'S'==temp[2]:
    w=3
   elif 'XX'==temp[2] or 'YY'==temp[2] or 'XY'==temp[2] or 'XZ'==temp[2]:
    w=3
   elif 'ZZ'==temp[2] or 'YZ'==temp[2]:
    w=3
   else:
    w=4

   n = len(temp) - w
   for z in range(0,n):
    #print z,index[z],x,w,temp
    j = int(index[z]) - 1
    #k = int(temp[0])
    #print x,j,w
    Smat[x,j] = temp[w]
    Smat[j,x] = temp[w]
    w = w + 1

if index_col2 != 0:
 line = negf_gamess2.readline()
 line2 = negf_gamess2.readline()
 line = negf_gamess2.readline()
 index = line2.split()
 h = len(index) 
 #print index
 si =numpy.int16(index[0])-1
 for x in range(si,natomic):
  line3 = negf_gamess2.readline()
  temp2 = line3.split()
  #print 't',x,temp
  if 'X'==temp2[2] or 'Y'==temp2[2] or 'Z'==temp2[2] or 'S'==temp2[2]:
   w=3
  elif 'XX'==temp2[2] or 'YY'==temp2[2] or 'XY'==temp2[2] or 'XZ'==temp2[2]:
   w=3
  elif 'ZZ'==temp2[2] or 'YZ'==temp2[2]:
   w=3
  else:
   w=4
  n = len(temp2) - w
  for z in range(0,n):
   j = int(index[z]) - 1
   k = int(temp2[0])
   Smat[x,j] = temp2[w]
   Smat[j,x] = temp2[w]
   w = w + 1


negf_smat = open(name + '_Smat','w')
for x in range(0,natomic):
 for y in range(0,natomic):
  #negf_smat.write("{0},".format(x) + "{0},".format(y) +  "{0}".format(Smat[x,y]) + "\n")
  negf_smat.write("{0}".format(Smat[x,y]) + "\n")


negf_smat.close()
negf_gamess2.close()
