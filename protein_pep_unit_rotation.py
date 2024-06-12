#!/usr/bin/python
from Ala_map_generator import Ala_MapsGenerator
from Gly_map_generator import Gly_MapsGenerator
import os
import sys
import smtplib
import datetime
from DbOperations import get_results_from_db
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders

def checkMultipleModel(pdb_id):
    infile=open("pdb_files/"+pdb_id+".pdb","r")
    towrite=""
    flag=0
    for line in infile:
        if line[0:5]=="MODEL":
            flag=1
            continue
        if line[0:4]=="ATOM" and flag==0:
            break
        if flag==1 and line[0:6]=="ENDMDL":
            break
        elif flag==1:
            towrite+=line
    infile.close()
    if flag!=0:
        opfile=open("pdb_files/"+pdb_id+".pdb","w")
        opfile.write(towrite)
        opfile.close()

def cleanPdbFile(pdb_id,chain_id,missing_res):
    infile=open("pdb_files/"+pdb_id+".pdb","r")
    cur_res=None
    residue=dict()
    towrite=[]
    for line in infile:
        if line[0:4]=="ATOM" and line[21]==chain_id:
            if cur_res==None:
                cur_res=int(line[22:26])
            elif cur_res!=int(line[22:26]) and cur_res not in missing_res: #new residue and not a missing residue, write N, CA and C atom lines, in that order (order important for phi-psi calculation)
                if residue.has_key("N"):
                    towrite.append(residue['N'])
                if residue.has_key("H"):
                    towrite.append(residue["H"])
                if residue.has_key("CD"):
                    towrite.append(residue["CD"])
                if residue.has_key("CA"):
                    towrite.append(residue['CA'])
                if residue.has_key("C"):
                    towrite.append(residue['C'])
                for item in residue.keys():
                    if item not in ['N','CA','C','H','CD']: #write rest of the atom lines of the residue
                        towrite.append(residue[item])
                residue=dict()
                cur_res=int(line[22:26])
            if cur_res in missing_res: #if missing residue, ignore
                residue=dict()
                cur_res=int(line[22:26])
            if cur_res==int(line[22:26]):
                if line[26]!=" ":
                    residue=dict()
                    cur_res=int(line[22:26])
                else:
                    if not residue.has_key(line[12:16].strip()):
                        residue[line[12:16].strip()]=line #store atom lines of non-missing residues
        elif line[0:6]=="ANISOU": #ignore anisou lines
            continue
        else:
            towrite.append(line)
    if residue.has_key("N"):
        towrite.append(residue['N'])
    if residue.has_key("H"):
        towrite.append(residue["H"])
    if residue.has_key("CD"):
        towrite.append(residue["CD"])
    if residue.has_key("CA"):
        towrite.append(residue['CA'])
    if residue.has_key("C"):
        towrite.append(residue['C'])
    for item in residue.keys():
        if item not in ['N','CA','C','H','CD']: #write rest of the atom lines of the residue
            towrite.append(residue[item])
    infile.close()
    opfile=open("pdb_files/"+pdb_id+".pdb","w") #rewrite existing pdb file
    for line in towrite:
        opfile.write(line)
    opfile.close()

def findMissingRes(pdb_id,chain_id):            
    infile=open(pdb_id+".pdb","r")
    cur_res=None
    missing_res=[]
    atoms=dict()
    for line in infile:
        if line[0:4] == "ATOM" and line[21]==chain_id:
            if cur_res==None:
                cur_res=int(line[22:26])
            if cur_res!=int(line[22:26]): #new residue
                if cur_res+1!=int(line[22:26]): #residue number currently read is not +1 of previously read residue
                    added_atoms=atoms.keys()
                    if "N" not in added_atoms or "C" not in added_atoms or "CA" not in added_atoms: #if any of the three backbone atoms not present in a residue, consider it as missing
                        missing_res.append(cur_res)
                    for i in range(cur_res+1,int(line[22:26])):
                        missing_res.append(i) 
                else:
                    added_atoms=atoms.keys()
                    if "N" not in added_atoms or "C" not in added_atoms or "CA" not in added_atoms: #if any of the three backbone atoms not present in a residue, consider it as missing
                        missing_res.append(cur_res)
                cur_res=int(line[22:26])
                atoms=dict()
            if line[21]==chain_id and line[12:16].strip() in ["N","CA","C"]:
                atom_type=line[12:16].strip()
                x=float(line[30:38])
                y=float(line[38:46])
                z=float(line[46:54])
                if not atoms.has_key(atom_type): #store first conformer (in case of multiple occupancy)
                    atoms[atom_type]=[x,y,z]
    infile.close()
    return missing_res

def ala_rotate(atoms_to_rotate,resname,resnumber,prepro):
    phipsi=Ala_MapsGenerator(atoms_to_rotate,pdb_id,chain_id,resname,resnumber,prepro)
    htmlfile.write('<tr>\n')
    htmlfile.write('<td>'+str(resnumber)+'</td>\n')
    htmlfile.write('<td>'+resname+'</td>\n')
    htmlfile.write('<td><a href="http://pauling.mbu.iisc.ac.in/bond_param_rmap/results/'+pdb_id+'_'+chain_id+'_'+str(resnumber)+'_rmap.png">Image file</a></td>\n')
    htmlfile.write('<td><a href="http://pauling.mbu.iisc.ac.in/bond_param_rmap/results/'+pdb_id+'_'+chain_id+'_'+resname+str(resnumber)+'.csv">Data file</a></td>\n')
    htmlfile.write('<td>'+str(phipsi[0])+'</td>\n')
    htmlfile.write('<td>'+str(phipsi[1])+'</td>\n')
    if phipsi[2]==1:
        htmlfile.write('<td>This residue\'s &straightphi;-&psi; is fully allowed according to its Ramachandran Map</td>\n')
    elif phipsi[2]==2:
        htmlfile.write('<td>This residue\'s &straightphi;-&psi; is partially allowed according to its Ramachandran Map</td>\n')
    else:
        htmlfile.write('<td>This residue\'s &straightphi;-&psi; is disallowed according to its Ramachandran Map</td>\n')
    htmlfile.write('</tr>\n')
def gly_rotate(atoms_to_rotate,resname,resnumber,prepro):
    phipsi=Gly_MapsGenerator(atoms_to_rotate,pdb_id,chain_id,resname,resnumber,prepro)
    htmlfile.write('<tr>\n')
    htmlfile.write('<td>'+str(resnumber)+'</td>\n')
    htmlfile.write('<td>'+resname+'</td>\n')
    htmlfile.write('<td><a href="http://pauling.mbu.iisc.ac.in/bond_param_rmap/results/'+pdb_id+'_'+chain_id+'_'+str(resnumber)+'_rmap.png">Image file</a></td>\n')
    htmlfile.write('<td><a href="http://pauling.mbu.iisc.ac.in/bond_param_rmap/results/'+pdb_id+'_'+chain_id+'_'+resname+str(resnumber)+'.csv">Data file</a></td>\n')
    htmlfile.write('<td>'+str(phipsi[0])+'</td>\n')
    htmlfile.write('<td>'+str(phipsi[1])+'</td>\n')
    if phipsi[2]==1:
        htmlfile.write('<td>This residue\'s &straightphi;-&psi; is fully allowed according to its Ramachandran Map</td>\n')
    elif phipsi[2]==2:
        htmlfile.write('<td>This residue\'s &straightphi;-&psi; is partially allowed according to its Ramachandran Map</td>\n')
    else:
        htmlfile.write('<td>This residue\'s &straightphi;-&psi; falls is disallowed according to its Ramachandran Map</td>\n')
    htmlfile.write('</tr>\n')
def send_to_queue(atom):
    global peptide_queue
    global atleast_one_map
    peptide_queue.append(atom)
    if atom[0]=="CA":
        atom_names=[row[0] for row in peptide_queue]
        ca_count=atom_names.count("CA")
        if ca_count==3:
            first_res=peptide_queue[0][5]
            last_res=peptide_queue[-1][5]
            if last_res-first_res!=2:
                new_pep=[]
                for item in peptide_queue:
                    if item[5]==last_res-1:
                        new_pep.append(item)
                    elif item[5]==last_res:
                        if item[0] in ["CA","C","O"]:
                            new_pep.append(item)
                peptide_queue=[]
                peptide_queue = new_pep[:]
            elif last_res-first_res==2:
                first_res_res=[]
                second_res_res=[]
                third_res_res=[]
                for item in peptide_queue:
                    if item[5]==first_res:
                        first_res_res.append(item[0])
                    elif item[5]==first_res+1:
                        second_res_res.append(item[0])
                        second_res_type=item[4]
                    else:
                        third_res_res.append(item[0])
                first_res_flag=second_res_flag=third_res_flag=0
                if "CA" in first_res_res and "C" in first_res_res and "O" in first_res_res:
                    first_res_flag=1
                if "N" in second_res_res and "H" in second_res_res and "CA" in second_res_res and "C" in second_res_res and "O" in second_res_res:
                    if second_res_type!="GLY":
                        if "CB" in second_res_res and "HA" in second_res_res:
                            second_res_flag=1
                    else:
                        if "HA2" in second_res_res and "HA3" in second_res_res:
                            second_res_flag=2
                if "N" in third_res_res and  "CA" in third_res_res:
                    if "H" in third_res_res:
                        third_res_flag=1
                    elif "CD" in third_res_res:
                        third_res_flag=2
                if second_res_flag==0:
                    if "CA" in second_res_res and "C" in second_res_res and "O" in second_res_res:
                        second_res_flag=-1
                if first_res_flag==1 and (second_res_flag==1 or second_res_flag==2) and (third_res_flag==1 or third_res_flag==2):
                    atoms_to_rotate=dict()
                    new_pep=[]
                    for item in peptide_queue:
                        if item[5]==first_res:
                            if item[0]=="CA":
                                atoms_to_rotate["CA1"]=[item[1],item[2],item[3]]
                            elif item[0]=="C":
                                atoms_to_rotate["C1"]=[item[1],item[2],item[3]]
                            elif item[0]=="O":
                                atoms_to_rotate["O1"]=[item[1],item[2],item[3]]
                        elif item[5]==first_res+1:
                            if item[0]=="N":
                                atoms_to_rotate["N2"]=[item[1],item[2],item[3]]
                            elif item[0]=="H":
                                atoms_to_rotate["H2"]=[item[1],item[2],item[3]]
                            elif item[0]=="CA":
                                atoms_to_rotate["CA2"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="CB":
                                atoms_to_rotate["CB"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="HA":
                                atoms_to_rotate["HA"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="HA2":
                                atoms_to_rotate["HA1"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="HA3":
                                atoms_to_rotate["HA2"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="C":
                                atoms_to_rotate["C2"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="O":
                                atoms_to_rotate["O2"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                        elif item[5]==last_res:
                            if item[0]=="N":
                                atoms_to_rotate["N3"]=[item[1],item[2],item[3]]
                            elif item[0]=="H":
                                atoms_to_rotate["H3"]=[item[1],item[2],item[3]]
                            elif item[0]=="CA":
                                atoms_to_rotate["CA3"]=[item[1],item[2],item[3]]
                            elif item[0]=="CD":
                                atoms_to_rotate["CD3"]=[item[1],item[2],item[3]]
                            new_pep.append(item)
                    peptide_queue=[]
                    peptide_queue=new_pep[:]
                    if second_res_flag==1:
                        if third_res_flag==1:
                            ala_rotate(atoms_to_rotate,second_res_type,first_res+1,"normal")
                            atleast_one_map=True
                        elif third_res_flag==2:
                            ala_rotate(atoms_to_rotate,second_res_type,first_res+1,"prepro")
                            atleast_one_map=True
                    else:
                        if third_res_flag==1:
                            gly_rotate(atoms_to_rotate,"GLY",first_res+1,"normal")
                            atleast_one_map=True
                        elif third_res_flag==2:
                            gly_rotate(atoms_to_rotate,"GLY",first_res+1,"prepro")
                            atleast_one_map=True
                elif second_res_flag==-1 and third_res_flag==1:
                    new_pep=[]
                    for item in peptide_queue:
                        if item[5]==first_res+1 and item[0] in ["CA","C","O","CB","HA"]:
                            new_pep.append(item)
                        elif item[5]==last_res and item[0] in ["CA","N","H"]:
                            new_pep.append(item)
                    peptide_queue=[]
                    peptide_queue=new_pep[:]
                elif first_res_flag==0 and second_res_flag==1 and third_res_flag==1:
                    new_pep=[]
                    for item in peptide_queue:
                        if item[5]==first_res+1 or item[5]==last_res:
                            new_pep.append(item)
                    peptide_queue=[]
                    peptide_queue=new_pep[:]
                else:
                    peptide_queue=[]
                    peptide_queue.append(atom)
                            
                    
           
                        
                        
logfile=open("rotation_log.log","w")                   
pdb_id=sys.argv[1]
chain_id=sys.argv[2]
present_in_db=sys.argv[3]
mail_id=sys.argv[4]
logfile.write("in pep unit program\n")
logfile.close()
if present_in_db=="True":
    logfile=open("rotation_log.log","a") 
    logfile.write("in presen in db if condition\n")
    logfile.close()
    uniq_id=get_results_from_db(pdb_id,chain_id)
else:
    missing_res=[]
    checkMultipleModel(pdb_id)
    logfile=open("rotation_log.log","a") 
    logfile.write("multiple models checked\n")
    logfile.close()
    cleanPdbFile(pdb_id,chain_id,missing_res)
    logfile=open("rotation_log.log","a") 
    logfile.write("file cleaned\n")
    logfile.close()
    first_n=True
    first_h=True
    global atleast_one_map
    atleast_one_map=False
    global peptide_queue
    peptide_queue=[]
    cur_date=datetime.datetime.now()
    uniq_id=str(cur_date.year)+str(cur_date.month)+str(cur_date.day)+"_"+str(cur_date.hour)+str(cur_date.minute)+str(cur_date.second)
    htmlfile=open('results/'+pdb_id+'_'+chain_id+'_'+uniq_id+'.html','w')
    htmlfile.write('<html>\n')
    htmlfile.write('<head>\n')
    htmlfile.write('<title>Bond-parameter specific Ramchandran Maps</title>\n')
    htmlfile.write('<link rel="stylesheet" type = "text/css" href = "parama.css" media = " all" />\n')
    htmlfile.write('<style>\n')
    htmlfile.write('td { border: 1px solid #ddd; padding:8px; background-color: #f2f2f2; text-align:center;}\n')
    htmlfile.write('th {background-color: #3399ff;border: 1px solid black;padding:8px;text-align:center;}\n')
    htmlfile.write('tr:nth-child(even) {background-color: #f2f2f2;}\n')
    htmlfile.write('</style>\n')
    htmlfile.write('</head>\n')
    htmlfile.write('<body>\n')
    htmlfile.write('<div class="header">\n')
    htmlfile.write('<h2 style="color:#3399ff;text-align:center;font-size:30px;">BORAM</h2>\n')
    htmlfile.write('<p style="font-size:25px;font-family: Verdana,Geneva, sans-serif; text-align:center;color:#3399ff;"><b style="color:#3399ff;">Bo</b>nd-parameter specific <b style="color:#3399ff;">Ra</b>machandran <b style="color:#3399ff;">M</b>aps </p>\n')
    htmlfile.write('</div>\n')
    htmlfile.write('<h3>Links to the generated Ramachandran maps and the corresponding data file are given in the following table</h3>\n')
    htmlfile.write('<table style="margin-left:auto;margin-right:auto;border-collapse:collapse;width:90%;">\n')
    htmlfile.write('<tr>\n')
    htmlfile.write('<th>Residue Number</th>\n')
    htmlfile.write('<th>Residue Name</th>\n')
    htmlfile.write('<th>Map Image file</th>\n')
    htmlfile.write('<th>Map Data</th>\n')
    htmlfile.write('<th>&straightphi; (&deg;)</th>\n')
    htmlfile.write('<th>&psi; (&deg;)</th>\n')
    htmlfile.write('<th>Fully Allowed/Partially Allowed/Disallowed</th>\n')
    htmlfile.write('</tr>\n')
    pdbfile=open("pdb_files/"+pdb_id+".pdb","r")
    logfile=open("rotation_log.log","a")
    logfile.write("Opened pdb file\n")
    logfile.close()
    for line in pdbfile:
        if line.startswith("ATOM"):
            if chain_id in line[21]:
                if first_n and line[12:16].strip()=="N":
                    first_n=False
                elif first_h and line[12:16].strip() in ["H1","H2","H3","H"]:
                    first_h=False
                elif line[12:16].strip()=="N":
                    atom_to_consider=["N",float(line[30:38]),float(line[38:46]),float(line[46:54]),line[17:20],int(line[22:26])]
                    send_to_queue(atom_to_consider)
                elif line[12:16].strip() in ["CA","C","O","CB","H","HA","HA2","HA3"]:
                    atom_to_consider=[line[12:16].strip(),float(line[30:38]),float(line[38:46]),float(line[46:54]),line[17:20],int(line[22:26])]
                    send_to_queue(atom_to_consider)
                if line[17:20]=="PRO" and line[12:16].strip()=="CD":
                    atom_to_consider=[line[12:16].strip(),float(line[30:38]),float(line[38:46]),float(line[46:54]),line[17:20],int(line[22:26])]
                    send_to_queue(atom_to_consider)
    pdbfile.close()
    htmlfile.write('</table>\n')
    htmlfile.write('</body>\n')
    htmlfile.write('</html>\n')
    htmlfile.close()
    if not atleast_one_map:
        cur_date=datetime.datetime.now()
        uniq_id=str(cur_date.year)+str(cur_date.month)+str(cur_date.day)+"_"+str(cur_date.hour)+str(cur_date.minute)+str(cur_date.second)
        htmlfile=open('results/'+pdb_id+'_'+chain_id+'_'+uniq_id+'.html','w')
        htmlfile.write('<html>\n')
        htmlfile.write('<head>\n')
        htmlfile.write('<title>Bond-parameter specific Ramchandran Maps</title>\n')
        htmlfile.write('<link rel="stylesheet" type = "text/css" href = "parama.css" media = " all" />\n')
        htmlfile.write('</head>\n')
        htmlfile.write('<body>\n')
        htmlfile.write('<div class="header">\n')
        htmlfile.write('<h2 style="color:#3399ff;text-align:center;font-size:30px;">BORAM</h2>\n')
        htmlfile.write('<p style="font-size:25px;font-family: Verdana,Geneva, sans-serif; text-align:center;color:#3399ff;"><b style="color:#3399ff;">Bo</b>nd-parameter specific <b style="color:#3399ff;">Ra</b>machandran <b style="color:#3399ff;">M</b>aps </p>\n')
        htmlfile.write('</div>\n')
        htmlfile.write('<h3> There were no two-linked peptide units with all the backbone atoms (including hydrogens) present.</h3>')
        htmlfile.close()
msg=MIMEMultipart()
toaddr=sys.argv[4]
msg['Subject']='Your results for the job in BoRam'
msg['From']='boram.server@gmail.com'
msg['To']=toaddr 
body="Dear User,\n\nPlease find attached the detailed results as well as the Ramachandran maps for the protein "+pdb_id+"_"+chain_id+" which you had submitted to BoRam. "
msg.attach(MIMEText(body,'plain'))
filename=pdb_id+"_"+chain_id+"_"+uniq_id+".html"
attachment=open("results/"+filename,"rb")
p=MIMEBase('application', 'octet-stream')
p.set_payload((attachment).read())
encoders.encode_base64(p)
p.add_header('Content-Disposition', "attachment; filename= %s" % filename)
msg.attach(p)  
server=smtplib.SMTP('smtp.gmail.com',587)
server.starttls()
server.login('boram.server@gmail.com','nsgrp2837')    
text=msg.as_string()
server.sendmail('boram.server@gmail.com',toaddr,text)
server.quit()       
