import math
import os

def make_sigma(sigma):
	new_sigma=10*(2**(1/6))*float(sigma)
	return new_sigma
def make_epsilon(epsilon):
	new_epsilon=float(epsilon)/-4.184
	return new_epsilon
def make_kbond(kbond):
	new_kbond=float(kbond)/836.8
	return new_kbond
def make_bonddistance(bdist):
	new_dist=float(bdist)*10
	return new_dist
def make_kangle(kang):
	new_kang=float(kang)/8.368
	return new_kang
def make_kdih(kdih):
	new_kdih=float(kdih)/4.184
	return new_kdih
def make_cmap(cmap):
	new_cmap=float(cmap)/4.184
	return new_cmap
SCAL_FACT=2

print("\n \nGMX2ORCA Version 1.0  \n \nThis program writes a generic GROMACS topology for any given forcefield in the ORCA QM/MM input format\n  \n The code requires an input.txt file that contains info about:\n 1) The structure of the system in .pdb format \n 2) The .itp files needed\n 3) How many molecules of a specific .itp are present \n 4) The number of atoms of the molecule \n 5) 6) 7) If the molecule has bonds,angles and dihedrals \n 8) If the .itp files need to be written with the ffbonded parameters\n\n\n-------------input file example-----------------\n\n;structure file name\nSystem.pdb\n;itp_name    nmols  natoms  dih(y/n)   angles(y/n)   bonds(y/n)  write_itp_prms(y/n)  cmaps(y/n)\nprotein.itp      1  5914        y          y           y                n                  y\nwater.itp    10000     3        n          y           y                n                  n\nions.itp        20     1        n          n           n                n                  n\n\n------------------------------------------------ \n \n \n NOTE 1: For the program to work correctly, the .pdb file needs to have the atom element as last column \n NOTE 2: For the program to work correctly, the .itp files should be in the same directory of the program, along with ffbonded.itp and ffnonbonded.itp \n NOTE 3: The program overwrites the .itp file that need to be manipulated as specified in the input.txt \n NOTE 4: The input .itp should be continous and the dihedral and improper dihedral section should be joined \n NOTE 5: If the program encounters an '#' either stops writing or encounters an error, make sure that there are no # in between section, The ';' is ignored by the program\n\n\n The program is running...")
file_mod=open("input.txt",'r').read()
righe_mod=file_mod.splitlines()
rm=3
lemon=0
ccmap=0
for riga in righe_mod:
	if riga=="":
		continue
	if riga.startswith(";"):
		continue
	elem=riga.split()
	if len(elem)<5:
		continue
	if elem[7]=="y":
		ccmap=ccmap+1
while rm<len(righe_mod):
	if righe_mod[rm]=="" or righe_mod[rm]==" ":
		continue
	elem=righe_mod[rm].split()
	if elem[6]=="n":
		rm=rm+1
		continue
	chcmap=elem[7]
	lemon=lemon+1
	fbonded=open("ffbonded.itp",'r').read()
	lines_fbonded=fbonded.splitlines()
	nin=elem[0]
	finput=open(nin,'r').read()
	itp_lines=finput.splitlines()
	index_atype=[]
	i=0
	atype=set()
	aname=set()
	while i<len(itp_lines):
		line=itp_lines[i]
		if line=="[ atoms ]":
			while True:
				line=itp_lines[i]
				elem=line.split()
				if line=="[ bonds ]":
					break
				if line.startswith(";"):
					i=i+1
					continue
				elif len(elem)<5:
					i=i+1
					continue
				else:
					index_atype.append([elem[0],elem[1]])
					atype.add(elem[1])
					aname.add(elem[4])
					i=i+1
		i=i+1
	i=0
	bonds=[]
	while i<len(itp_lines):
		line=itp_lines[i]
		if line=="[ bonds ]":
			i=i+1
			while True:
				line=itp_lines[i]
				elem=line.split()
				if line=="[ angles ]" or line=="[ exclusions ]" or line=="[ pairs ]":
					break
				if line.startswith(";"):
					i=i+1
					continue
				if len(line)<1:
					i=i+1
					continue
				else:
					bonds.append([elem[0],elem[1]])
					i=i+1
		i=i+1
	i=0
	ffbonded_bonds=[]
	while i<len(bonds):
		nesimo_bond=bonds[i]
		conv_bond=[]
		j=0
		counter=0
		while counter<2:
			if j>=len(index_atype):
				j=0
				continue
			coppia=index_atype[j]
			if nesimo_bond[0]==coppia[0]:
				conv_bond.append(coppia[1])
				counter=counter+1  
			elif counter==1 and nesimo_bond[1]==coppia[0]:
				conv_bond.append(coppia[1])
				counter=counter+1
			j=j+1
		ffbonded_bonds.append(conv_bond)
		i=i+1
		
	i=0
	ffbonded_truebonds=[]
	while i<len(lines_fbonded):
		line=lines_fbonded[i]
		if line=="[ bondtypes ]":
			i=i+1
			while True:
				line=lines_fbonded[i]
				elem=line.split()
				if line=="[ constrainttypes ]" or line.startswith("#") or line=="[ angletypes ]":
					break
				if line.startswith(";"):
					i=i+1
					continue
				elif len(elem)<2:
					i=i+1
					continue
				else:
					if (elem[0] in atype) and (elem[1] in atype):
						ffbonded_truebonds.append(line)
					i=i+1
		i=i+1	
	i=0
	bonds_correct=[]
	while i<len(bonds):
		ind1=bonds[i][0]
		ind2=bonds[i][1]
		t1=ffbonded_bonds[i][0]
		t2=ffbonded_bonds[i][1]
		j=0
		while j<len(ffbonded_truebonds):
			elem=ffbonded_truebonds[j].split()
			if (t1==elem[0] or t1==elem[1]) and (t2==elem[0] or t2==elem[1]):
				funct=elem[2]
				dist=elem[3]
				kbond=elem[4]
				bonds_correct.append([ind1,ind2,funct,dist,kbond])
				break
			j=j+1
		i=i+1

	fout=open(nin,'w')
	change=False
	k=0
	for line in itp_lines:
		if k>=len(itp_lines):
			break
		if line=="[ bonds ]":
			fout.write(line+"\n")
			fout.write(itp_lines[k+1]+"\n")
			j=0
			k=k+2
			while j<len(bonds_correct):
				elem=bonds_correct[j]
				fout.write("{:7d}{:7d}{:3d}    {:<13.6f}{:<11.6f}".format(int(elem[0]),int(elem[1]),int(elem[2]),float(elem[3]),float(elem[4]))+"\n")
				j=j+1
				k=k+1
			change=True
		if change==False:
			fout.write(line+"\n")
			k=k+1
		else:
			fout.write(itp_lines[k]+"\n")
			k=k+1
	del(bonds)
	del(bonds_correct)
	del(ffbonded_bonds)
	del(ffbonded_truebonds)
	fout.close()
	#################################################################################### RISCRITTA LA SEZIONE BONDS
	angles=[]
	i=0
	while i<len(itp_lines):
		line=itp_lines[i]
		if line=="[ angles ]":
			i=i+1
			while True:
				line=itp_lines[i]
				elem=line.split()
				if line=="[ dihedrals ]" or line=="[ exclusions ]" or line=="[ pairs ]":
					break
				if line.startswith(";"):
					i=i+1
					continue
				if len(line)<1:
					i=i+1
					continue
				else:
					angles.append([elem[0],elem[1],elem[2]])
					i=i+1
		i=i+1
	i=0
	ffbonded_angles=[]
	while i<len(angles):
		nesimo_angle=angles[i]
		conv_angle=[]
		j=0
		counter=0
		while counter<3:
			if j>=len(index_atype):
				j=0
				continue
			coppia=index_atype[j]
			if nesimo_angle[0]==coppia[0]:
				conv_angle.append(coppia[1])
				counter=counter+1  
			elif counter==1 and nesimo_angle[1]==coppia[0]:
				conv_angle.append(coppia[1])
				counter=counter+1
			elif counter==2 and nesimo_angle[2]==coppia[0]:
				conv_angle.append(coppia[1])
				counter=counter+1			
			j=j+1
		ffbonded_angles.append(conv_angle)
		i=i+1
		
	i=0
	ffbonded_trueangles=[]
	while i<len(lines_fbonded):
		line=lines_fbonded[i]
		if line=="[ angletypes ]":
			i=i+1
			while True:
				line=lines_fbonded[i]
				elem=line.split()
				if line=="[ constrainttypes ]" or line.startswith("#") or line=="[ angletypes ]" or line=="[ dihedraltypes ]":
					break
				if line.startswith(";"):
					i=i+1
					continue
				elif len(elem)<2:
					i=i+1
					continue
				else:
					if (elem[0] in atype) and (elem[1] in atype) and (elem[2] in atype):
						ffbonded_trueangles.append(line)
					i=i+1
		i=i+1	
	i=0
	angles_correct=[]
	while i<len(angles):
		ind1=angles[i][0]
		ind2=angles[i][1]
		ind3=angles[i][2]
		t1=ffbonded_angles[i][0]
		t2=ffbonded_angles[i][1]
		t3=ffbonded_angles[i][2]
		j=0
		while j<len(ffbonded_trueangles):
			elem=ffbonded_trueangles[j].split()
			if (t1==elem[0] or t1==elem[1] or t1==elem[2]) and (t2==elem[0] or t2==elem[1] or t2==elem[2]) and (t3==elem[0] or t3==elem[1] or t3==elem[2]):
				funct=elem[3]
				ang=elem[4]
				kang=elem[5]
				angles_correct.append([ind1,ind2,ind3,funct,ang,kang])
				break
			j=j+1
		i=i+1
	finput=open(nin,'r').read()
	itp_lines=finput.splitlines()
	fout=open(nin,'w')
	change=False
	k=0
	for line in itp_lines:
		if k>=len(itp_lines):
			break
		if line=="[ angles ]":
			fout.write(line+"\n")
			fout.write(itp_lines[k+1]+"\n")
			j=0
			k=k+2
			while j<len(angles_correct):
				elem=angles_correct[j]
				fout.write("{:7d}{:7d}{:7d}{:3d}    {:<13.6f}{:<11.6f}".format(int(elem[0]),int(elem[1]),int(elem[2]),int(elem[3]),float(elem[4]),float(elem[5]))+"\n")
				j=j+1
				k=k+1
			change=True
		if change==False:
			fout.write(line+"\n")
			k=k+1
		else:
			fout.write(itp_lines[k]+"\n")
			k=k+1
	del(angles)
	del(angles_correct)
	del(ffbonded_angles)
	del(ffbonded_trueangles)
	############################################################################################## RISCRITTI GLI ANGLES
	dihedrals=[]
	i=0
	while i<len(itp_lines):
		line=itp_lines[i]
		if line=="[ dihedrals ]":
			i=i+1
			while True:
				line=itp_lines[i]
				elem=line.split()
				if line=="[ dihedrals ]" or line=="[ exclusions ]" or line=="[ pairs ]" or line.startswith("#") or line=="[ cmap ]" or line=="" :
					break
				if line.startswith(";"):
					i=i+1
					continue
				if len(line)<1:
					i=i+1
					continue
				else:
					dihedrals.append([elem[0],elem[1],elem[2],elem[3]])
					i=i+1
		i=i+1
		
	i=0
	ffbonded_dihedrals=[]
	while i<len(dihedrals):
		nesimo_dihedral=dihedrals[i]
		conv_dihedral=[]
		j=0
		counter=0
		while counter<4:
			if j>=len(index_atype):
				j=0
				continue
			coppia=index_atype[j]
			if counter==0 and nesimo_dihedral[0]==coppia[0]:
				conv_dihedral.append(coppia[1])
				counter=counter+1  
			elif counter==1 and nesimo_dihedral[1]==coppia[0]:
				conv_dihedral.append(coppia[1])
				counter=counter+1
			elif counter==2 and nesimo_dihedral[2]==coppia[0]:
				conv_dihedral.append(coppia[1])
				counter=counter+1
			elif counter==3 and nesimo_dihedral[3]==coppia[0]:
				conv_dihedral.append(coppia[1])
				counter=counter+1						
			j=j+1
		ffbonded_dihedrals.append(conv_dihedral)
		i=i+1
		
	i=0
	ffbonded_truedihedrals=[]
	mask="X"
	while i<len(lines_fbonded):
		line=lines_fbonded[i]
		if line=="[ dihedraltypes ]":
			i=i+1
			while True:
				line=lines_fbonded[i]
				elem=line.split()
				if line=="[ constrainttypes ]" or line.startswith("#") or line=="[ dihedraltypes ]" or line=="[ dihedraltypes ]" or line=="" or line=="[ cmap ]":
					break
				if line.startswith(";"):
					i=i+1
					continue
				elif len(elem)<2:
					i=i+1
					continue
				else:
					if len(elem)==7:
						line=line+"   1"
					if (elem[0] in atype) and (elem[1] in atype) and (elem[2] in atype) and (elem[3] in atype):
						ffbonded_truedihedrals.append(line)
					elif elem[0]==mask or elem[1]==mask or elem[2]==mask or elem[3]==mask:
						ts=[elem[0],elem[1],elem[2],elem[3]]
						if any(tx in atype and ty in atype for tx in ts for ty in ts if tx !=ty):
							ffbonded_truedihedrals.append(line)			
					i=i+1
		i=i+1	
	i=0
	dihedrals_correct=[]
	while i<len(dihedrals):
		ind1=dihedrals[i][0]
		ind2=dihedrals[i][1]
		ind3=dihedrals[i][2]
		ind4=dihedrals[i][3]
		t1=ffbonded_dihedrals[i][0]
		t2=ffbonded_dihedrals[i][1]
		t3=ffbonded_dihedrals[i][2]
		t4=ffbonded_dihedrals[i][3]
		j=0
		while j<len(ffbonded_truedihedrals):
			elem=ffbonded_truedihedrals[j].split()
			if (t1==elem[0] or t1==elem[1] or t1==elem[2] or t1==elem[3]) and (t2==elem[0] or t2==elem[1] or t2==elem[2] or t2==elem[3]) and (t3==elem[0] or t3==elem[1] or t3==elem[2] or t3==elem[3]) and (t4==elem[0] or t4==elem[1] or t4==elem[2] or t4==elem[3]):
				funct=elem[4]
				dih=elem[5]
				kdih=elem[6]
				period=elem[7]
				dihedrals_correct.append([ind1,ind2,ind3,ind4,funct,dih,kdih,period])
				break
			elif elem[0]==mask or elem[1]==mask or elem[2]==mask or elem[3]==mask:
				ts=[t1,t2,t3,t4]
				if any(tx in atype and ty in atype for tx in ts for ty in ts if tx !=ty):
					funct=elem[4]
					dih=elem[5]
					kdih=elem[6]
					period=elem[7]			
					dihedrals_correct.append([ind1,ind2,ind3,ind4,funct,dih,kdih,period])	
					break		
			j=j+1
		i=i+1
	fout.close()
	finput=open(nin,'r').read()
	itp_lines=finput.splitlines()
	fout=open(nin,'w')
	change=False
	k=0
	for line in itp_lines:
		if k>=len(itp_lines):
			break
		if line=="[ dihedrals ]":
			fout.write(line+"\n")
			fout.write(itp_lines[k+1]+"\n")
			j=0
			k=k+2
			while j<len(dihedrals_correct):
				elem=dihedrals_correct[j]
				fout.write("{:7d}{:7d}{:7d}{:7d}{:3d}    {:<13.6f}{:<11.6f} {:7d}".format(int(elem[0]),int(elem[1]),int(elem[2]),int(elem[3]),int(elem[4]),float(elem[5]),float(elem[6]),int(elem[7]))+"\n")
				j=j+1
				k=k+1
			change=True
		if change==False:
			fout.write(line+"\n")
			k=k+1
		else:	
			if k>=len(itp_lines):
				break
			fout.write(itp_lines[k]+"\n")
			k=k+1
	del(dihedrals)
	del(dihedrals_correct)
	del(ffbonded_dihedrals)
	del(ffbonded_truedihedrals)
	fout.close()
#################################################################################################################### RISCRITTI I DIHEDRALS
	cmaps=[]
	num_cmaps=0
	if chcmap=="n":
		rm=rm+1
		continue
	i=0
	while i<len(itp_lines):
		line=itp_lines[i]
		if line=="[ cmap ]":
			i=i+1
			while True:
				line=itp_lines[i]
				elem=line.split()
				if line=="[ dihedrals ]" or line=="[ exclusions ]" or line=="[ pairs ]" or line=="" or line.startswith("#"):
					break
				if line.startswith(";"):
					i=i+1
					continue
				if len(line)<1:
					i=i+1
					continue
				else:
					cmaps.append([elem[0],elem[1],elem[2],elem[3],elem[4]])
					num_cmaps=num_cmaps+1
					i=i+1
		i=i+1
	i=0
	ffbonded_cmaps=[]
	while i<len(cmaps):
		nesimo_cmap=cmaps[i]
		conv_cmap=[]
		j=0
		counter=0
		while counter<5:
			if j>=len(index_atype):
				j=0
				continue
			coppia=index_atype[j]
			if nesimo_cmap[0]==coppia[0]:
				conv_cmap.append(coppia[1])
				counter=counter+1  
			elif counter==1 and nesimo_cmap[1]==coppia[0]:
				conv_cmap.append(coppia[1])
				counter=counter+1
			elif counter==2 and nesimo_cmap[2]==coppia[0]:
				conv_cmap.append(coppia[1])
				counter=counter+1
			elif counter==3 and nesimo_cmap[2]==coppia[0]:
				conv_cmap.append(coppia[1])
				counter=counter+1		
			elif counter==4 and nesimo_cmap[2]==coppia[0]:
				conv_cmap.append(coppia[1])
				counter=counter+1		
			j=j+1
		ffbonded_cmaps.append(conv_cmap)
		i=i+1
	i=0
	fcmap=open("cmap.itp",'r').read()
	lines_fbonded=fcmap.splitlines()
	lines_fbonded.append("###")
	ffbonded_truecmaps=[]
	cmaps_val=[]
	while i<len(lines_fbonded):
		line=lines_fbonded[i]
		if line=="[ cmaptypes ]":
			i=i+1
			while True:
				line=lines_fbonded[i]
				elem=line.split()
				if line=="":
					i=i+1
					continue
				if line=="[ constrainttypes ]" or line=="###" or line=="[ cmaptypes ]" or line=="[ dihedraltypes ]":
					break
				if line.startswith(";"):
					i=i+1
					continue
				elif len(elem)<2:
					i=i+1
					continue
				else:
					if (elem[0] in atype) and (elem[1] in atype) and (elem[2] in atype) and (elem[3] in atype) and (elem[4] in atype):
						ffbonded_truecmaps.append([elem[0],elem[1],elem[2],elem[3],elem[4]])
						r=0
						i=i+1
						single_cmap=[]
						while r<58:
							line=lines_fbonded[i][:-2]
							single_cmap.append(line)
							i=i+1
							r=r+1
						cmaps_val.append(single_cmap)
						continue
					else:
						i=i+1
		i=i+1	
	i=0
	cmaps_correct=[]
	cmaps_val_correct=[]
	foutcmap=open("cmap.tmp",'w')
	while i<len(cmaps):
		ind1=cmaps[i][0]
		ind2=cmaps[i][1]
		ind3=cmaps[i][2]
		ind4=cmaps[i][3]
		ind5=cmaps[i][4]
		t1=ffbonded_cmaps[i][0]
		t2=ffbonded_cmaps[i][1]
		t3=ffbonded_cmaps[i][2]
		t4=ffbonded_cmaps[i][3]
		t5=ffbonded_cmaps[i][4]
		ts=[t1,t2,t3,t4,t5]
		j=0
		while j<len(ffbonded_truecmaps):
			elem=ffbonded_truecmaps[j]
			if all(tx in elem[:4] for tx in ts):
				foutcmap.write(str(ind1)+" "+str(ind2)+" "+str(ind3)+" "+str(ind4)+" "+str(ind5))
				#cmaps_correct.append([ind1,ind2,ind3,ind4,ind5])
				#cmaps_val_correct.append(cmaps_val[j])
				r=0
				stringona=" "
				cmapsval=cmaps_val[j]
				while r<len(cmapsval):
					val=cmapsval[r].split()
					for value in val:
						stringona=stringona+"   "+value
					r=r+1
				foutcmap.write(stringona+"\n")
				break
			j=j+1
		i=i+1
	foutcmap.close()
	del(cmaps_val)
	del(single_cmap)
	del(cmaps_val_correct)
	del(cmaps)
	del(ffbonded_cmaps)
	del(ffbonded_truecmaps)
	rm=rm+1
	fout.close()
if lemon != 0:
	print("Wrote itp files")
else:
	print("Skipping itp files manipulation")
#########################################################################################################################################
############################################################################ dal file ffbondedmod.py
#########################################################################################################################################
f_input=open("input.txt",'r').read()
lines_input=f_input.splitlines()
nstruc=lines_input[1]
alm_pdb=nstruc.split(".")
pdb_name=alm_pdb[0]
fin_struc=open(nstruc,'r')
fin_nonbonded=open("ffnonbonded.itp",'r').read()
righe_ffnonbonded=fin_nonbonded.splitlines()
elements=[]
anames=set()
atype=set()
natoms=0
fout2=open("names.txt",'w')
cc_max=0
for line in fin_struc:
	elem=line.split()
	if elem[0]=="ATOM" or elem[0]=="HETATM":
		natoms=natoms+1
		fout2.write(str(int(elem[1])+100000*cc_max)+" "+elem[2]+"\n")
		if int(elem[1])==99999:
			cc_max=cc_max+1
		elements.append(elem[-1])
		anames.add(elem[2])
fout=open("orca.prms",'w')
fout2.close()
fin_indexes=open("names.txt",'r').read()
righe_indexes=fin_indexes.splitlines()
i=0
indexes=[]
while i<len(righe_indexes):
	elem=righe_indexes[i].split()
	indexes.append(str(elem[0]))
	i=i+1
fout.write("$fftype"+"\n")
fout.write("CHARMM"+"\n")
fout.write("$atoms"+"\n")
fout.write(str(natoms)+" 1 4"+"\n")
f_input=open("input.txt",'r').read()
ind_itp=[]
############################################################ FINO A QUI ASSOCIA ATOMNAME AD ATOMTYPE
righe_input=f_input.splitlines()
r_in=3
aname_atype=[]
charges=[]
while r_in<len(righe_input):
	elem=righe_input[r_in].split()
	f_itp=open(elem[0],'r').read()
	righe_itp=f_itp.splitlines()
	nmols=int(elem[1])
	jj=0
	while jj<nmols:
		c=0
		j=0
		while j<len(righe_itp):
			if righe_itp[j]=="[ atoms ]":
				j=j+2
				while righe_itp[j]!="[ bonds ]":
					elem=righe_itp[j].split()
					if len(elem)<5:
						break
					if elem[0]==";":
						j=j+1
						continue
					if elem[4] in anames:
						aname_atype.append([elem[4],elem[1]])
						charges.append(elem[6])
						ind_itp.append(elem[0])
						atype.add(elem[1])
					else:
						c=c+1
					j=j+1
			j=j+1
		jj=jj+1
	r_in=r_in+1
i=-1
counter_index=0
for coppia in aname_atype:
	i=i+1
	for line in righe_ffnonbonded:
		elem=line.split()
		if len(elem)<1:
			continue
		if coppia[1]==elem[0]:
			sigma=elem[5]
			epsilon=elem[6]
			break
	rmin=make_sigma(sigma)
	epsilon=make_epsilon(epsilon)
	charge=charges[i]
	if float(charge)<0:
		fout.write("{:>6d}{:>4}{:>14.6f}{:>13.6f}{:>13.6f}{:>13.6f}{:>13.6f}".format(int(indexes[counter_index]),str(elements[i]),float(charge),epsilon,rmin,epsilon/SCAL_FACT,rmin)+"\n")
	else:
		fout.write("{:>6d}{:>4}{:>14.6f}{:>13.6f}{:>13.6f}{:>13.6f}{:>13.6f}".format(int(indexes[counter_index]),str(elements[i]),float(charge),float(epsilon),float(rmin),float(epsilon/SCAL_FACT),float(rmin))+"\n")
	counter_index=counter_index+1
############################################################### QUI SCRIVE LA SEZIONE DI CHARGES,SIGMA e EPSILON
print("Wrote Non-Bonded")

fout.write("$bonds"+"\n")
fout.write("nbonds 2 2"+"\n")
nbonds=0
fout2.close()
fin_indexes=open("names.txt",'r').read()
righe_indexes=fin_indexes.splitlines()
i=0
indexes=[]
while i<len(righe_indexes):
	elem=righe_indexes[i].split()
#	aname_atype[i].append(elem[0])
	indexes.append(elem[0])
	i=i+1
#aname_atype_index=aname_atype.copy()
#aname_atype=[]
counter_index=0
counter_bonds=0
f_input=open("input.txt",'r').read()
righe_input_itp=f_input.splitlines()
qxy=0
for file_itp in righe_input_itp:
	if qxy<3:
		qxy=qxy+1
		continue
	elem=file_itp.split()
	if elem[5]=='n':
		counter_index=counter_index+int(elem[2])
		continue
	f_itp=open(elem[0],'r').read()
	lines_itp=f_itp.splitlines()
	lines_itp = [line.replace('\t', ' ') for line in lines_itp]
	j=0
	while j<len(lines_itp):
		if lines_itp[j]=="[ bonds ]":
			j=j+2
			bonds=[]
			i_indexes=[]
			while True:
				if (lines_itp[j] == "[ angles ]") or (lines_itp[j] == "[ pairs ]") or (lines_itp[j] == "[ exclusions ]"):
					break
				check=True
				if len(lines_itp[j])>2:
					bonds.append(lines_itp[j])
					q=lines_itp[j].split()
					i_indexes.append(int(q[0]))
				j=j+1
			break
		j=j+1
	seen=set()
	i_indexes=[x for x in i_indexes if not (x in seen or seen.add(x))]
	i_indexes.sort()	
	nmols=int(elem[1])
	natoms_itp=int(elem[2])
	j=0
	while j<nmols:
		iitp=1
		kk=0
		while iitp<=natoms_itp:
			if iitp not in i_indexes:
				counter_index=counter_index+1
				iitp=iitp+1
				continue
			kkk=0
			while kkk<len(bonds):
				bond=bonds[kkk].split()
				if int(bond[0])==i_indexes[kk]:
					ind1=int(bond[0])
					ind2=int(bond[1])
					diff=ind2-ind1
					i_orca=indexes[counter_index]
					j_orca=indexes[counter_index+diff]
					kdist=make_bonddistance(bond[3])
					kbond=make_kbond(bond[4])
					fout.write("{:>6d}{:>11d}{:>11.6f} {:<11.6f}".format(int(i_orca),int(j_orca),float(kdist),float(kbond))+"\n")### DA FORMATTARE
					counter_bonds=counter_bonds+1
				kkk=kkk+1
			counter_index=counter_index+1
			kk=kk+1
			iitp=iitp+1
		j=j+1
print("Wrote Bonds")
################################################################QUI HA SCRITTO LA SEZIONE BONDS
fout.write("$angles"+"\n")
fout.write("nangles 3 2"+"\n")
counter_index=0
counter_angles=0
f_input=open("input.txt",'r').read()
righe_input_itp=f_input.splitlines()
qxy=0
for file_itp in righe_input_itp:
	if qxy<3:
		qxy=qxy+1
		continue
	elem=file_itp.split()
	if elem[4]=='n':
		counter_index=counter_index+int(elem[2])
		continue
	f_itp=open(elem[0],'r').read()
	lines_itp=f_itp.splitlines()
	lines_itp = [line.replace('\t', ' ') for line in lines_itp]
	j=0
	while j<len(lines_itp):
		if lines_itp[j]=="[ angles ]":
			j=j+2
			angles=[]
			i_indexes=[]
			while True:
				if (lines_itp[j] == "[ dihedrals ]") or (lines_itp[j] == "[ pairs ]") or (lines_itp[j] == "[ exclusions ]") or (lines_itp[j] == ""):
					break
				check=True
				if len(lines_itp[j])>2:
					angles.append(lines_itp[j])
					q=lines_itp[j].split()
					i_indexes.append(int(q[0]))
				j=j+1
			break
		j=j+1
	seen=set()
	i_indexes=[x for x in i_indexes if not (x in seen or seen.add(x))]
	i_indexes.sort()	
	nmols=int(elem[1])
	natoms_itp=int(elem[2])
	j=0
	while j<nmols:
		iitp=1
		kk=0
		while iitp<=natoms_itp:
			if iitp not in i_indexes:
				counter_index=counter_index+1
				iitp=iitp+1
				continue
			kkk=0
			while kkk<len(angles):
				angle=angles[kkk].split()
				if int(angle[0])==i_indexes[kk]:
					ind1=int(angle[0])
					ind2=int(angle[1])
					ind3=int(angle[2])
					diff1=ind2-ind1
					diff2=ind3-ind1
					i_orca=indexes[counter_index]
					j_orca=indexes[counter_index+diff1]
					k_orca=indexes[counter_index+diff2]
					ang=(angle[4])
					kang=make_kangle(angle[5])
					fout.write("{:>6d}{:>9d}{:>9d}{:>13.6f}{:>13.6f}".format(int(i_orca),int(j_orca),int(k_orca),float(ang),float(kang))+"\n")### DA FORMATTARE
					counter_angles=counter_angles+1
				kkk=kkk+1
			counter_index=counter_index+1
			kk=kk+1
			iitp=iitp+1
		j=j+1
print("Wrote Angles")
##################################################################QUI HA SCRITTO GLI ANGOLI
fout.write("$dihedrals"+"\n")
fout.write("ndih 4 3"+"\n")
counter_index=0
counter_dih=0
f_input=open("input.txt",'r').read()
righe_input_itp=f_input.splitlines()
qxy=0
for file_itp in righe_input_itp:
	if qxy<3:
		qxy=qxy+1
		continue
	elem=file_itp.split()
	if elem[3]=='n':
		counter_index=counter_index+int(elem[2])
		continue
	f_itp=open(elem[0],'r').read()
	lines_itp=f_itp.splitlines()
	lines_itp = [line.replace('\t', ' ') for line in lines_itp]
	j=0
	while j<len(lines_itp):
		if lines_itp[j]=="[ dihedrals ]":
			j=j+2
			dihedrals=[]
			i_indexes=[]
			while True:
				if (lines_itp[j] == "[ dihedrals ]") or (lines_itp[j] == "[ pairs ]") or (lines_itp[j] == "[ exclusions ]") or lines_itp[j] == "":
					break
				if len(lines_itp[j])>2:
					dihedrals.append(lines_itp[j])
					q=lines_itp[j].split()
					i_indexes.append(int(q[0]))
				j=j+1
			break
		j=j+1
	seen=set()
	i_indexes=[x for x in i_indexes if not (x in seen or seen.add(x))]
	i_indexes.sort()	
	nmols=int(elem[1])
	natoms_itp=int(elem[2])
	j=0
	while j<nmols:
		iitp=1
		kk=0
		while iitp<=natoms_itp:
			if iitp not in i_indexes:
				counter_index=counter_index+1
				iitp=iitp+1
				continue
			kkk=0
			while kkk<len(dihedrals):
				dihedral=dihedrals[kkk].split()
				if int(dihedral[0])==i_indexes[kk]:
					ind1=int(dihedral[0])
					ind2=int(dihedral[1])
					ind3=int(dihedral[2])
					ind4=int(dihedral[3])
					diff1=ind2-ind1
					diff2=ind3-ind1
					diff3=ind4-ind1
					i_orca=indexes[counter_index]
					j_orca=indexes[counter_index+diff1]
					k_orca=indexes[counter_index+diff2]
					l_orca=indexes[counter_index+diff3]
					dih=(dihedral[5])
					kdih=make_kdih(dihedral[6])
					periodicity=dihedral[7]
					fout.write("{:>6d}{:>9d}{:>9d}{:>9d}{:>13.6f}{:>13.6f}{:>9d}".format(int(i_orca),int(j_orca),int(k_orca),int(l_orca),float(dih),float(kdih),int(periodicity))+"\n")### DA FORMATTARE
					counter_dih=counter_dih+1
				kkk=kkk+1
			counter_index=counter_index+1
			kk=kk+1
			iitp=iitp+1
		j=j+1
print("Wrote Dihedrals")
################################################################################################################## SCRITTI I DIEDRI
fout.write("$cmap"+"\n")
fout.write("ncmap 8 576"+"\n")
fcmap=open("cmap.tmp",'r').read()
cmaps=fcmap.splitlines()
counter_index=0
counter_cmp=0
boundary=0
cccmap=0
qxy=0
mmm=0
i_indexes=[]
for file_itp in righe_input_itp:
	if file_itp=="":
		continue
	if qxy<3:
		qxy=qxy+1
		continue
	elem=file_itp.split()
	if elem[7]=='n':
		counter_index=counter_index+int(elem[2])
		continue
	f_itp=open(elem[0],'r').read()
	j=0
	numb_cmaps=0
	righe_itp=f_itp.splitlines()
	zzz=0
	for riga in righe_itp:
		if riga=="[ cmap ]":
			zzz=zzz+2
			while True:
				riga=righe_itp[zzz]
				if riga==" " or riga=="" or riga.startswith("#"):
					break
				numb_cmaps=numb_cmaps+1
				zzz=zzz+1
		zzz=zzz+1
	range_cmap=numb_cmaps+boundary
	subset_cmaps=[]
	while mmm<range_cmap:
		bbb=cmaps[mmm].split()
		subset_cmaps.append(cmaps[mmm])
		i_indexes.append(float(bbb[0]))
		mmm=mmm+1
		boundary=numb_cmaps
	cccmap=cccmap+1
	seen=set()
	i_indexes=[x for x in i_indexes if not (x in seen or seen.add(x))]
	i_indexes.sort()	
	nmols=int(elem[1])
	natoms_itp=int(elem[2])
	j=0
	while j<nmols:
		iitp=1
		kk=0
		while iitp<=natoms_itp:
			if iitp not in i_indexes:
				counter_index=counter_index+1
				iitp=iitp+1
				continue
			kkk=0
			while kkk<len(subset_cmaps):
				cmap=subset_cmaps[kkk].split()
				if int(cmap[0])==i_indexes[kk]:
					ind1=int(cmap[0])
					ind2=int(cmap[1])
					ind3=int(cmap[2])
					ind4=int(cmap[3])
					ind5=int(cmap[4])
					diff1=ind2-ind1
					diff2=ind3-ind1
					diff3=ind4-ind1
					diff4=ind5-ind1
					i_orca=indexes[counter_index]
					j_orca=indexes[counter_index+diff1]
					k_orca=indexes[counter_index+diff2]
					l_orca=indexes[counter_index+diff3]
					m_orca=indexes[counter_index+diff4]
					fout.write(" {} {} {} {} {} {} {} {}".format(int(i_orca),int(j_orca),int(k_orca),int(l_orca),(j_orca),int(k_orca),int(l_orca),int(m_orca))+"\n")### DA FORMATTARE
					counter_cmp=counter_cmp+1
					zzz=0
					r=0
					lmn=5
					while r<96:
						fout.write("{:>13.8f}{:>13.8f}{:>13.8f}{:>13.8f}{:>13.8f}{:>13.8f}".format(float(make_cmap(cmap[lmn])),float(make_cmap(cmap[lmn+1])),float(make_cmap(cmap[lmn+2])),float(make_cmap(cmap[lmn+3])),float(make_cmap(cmap[lmn+4])),float(make_cmap(cmap[lmn+5])))+"\n")
						lmn=lmn+6
						r=r+1
				kkk=kkk+1
			counter_index=counter_index+1
			kk=kk+1
			iitp=iitp+1
		j=j+1
fout.close()
finprov=open("orca.prms",'r').readlines()
fout=open(pdb_name+"_ORCAFF.prms",'w')
for line in finprov:
	elem=line.split()
	if len(elem)<1:
		fout.write(line)
		continue
	if elem[0]=="nbonds":
		new_line=str(counter_bonds)+line[len("nbonds"):]
		fout.write(new_line)
		continue
	elif elem[0]=="nangles":
		new_line=str(counter_angles)+line[len("nangles"):]
		fout.write(new_line)
		continue
	elif elem[0]=="ndih":
		new_line=str(counter_dih)+line[len("ndih"):]
		fout.write(new_line)
		continue
	elif elem[0]=="ncmap":
		new_line=str(counter_cmp)+line[len("ncmap"):]
		fout.write(new_line)
		continue
	fout.write(line)
os.remove("orca.prms")
#os.remove("cmap.tmp")
	
	
