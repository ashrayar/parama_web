# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 15:07:19 2018

@author: ashraya
"""
import MySQLdb
db = MySQLdb.connect("localhost","root","mysql123","boram" )
cursor = db.cursor()
#infile=open("chains_with_h.txt","r")
#for line in infile:
#    lineparts=line.split(",")
#    pdb_id=lineparts[0]
#    chain_id=lineparts[1]
#    rwork=float(lineparts[2])
#    rfree=float(lineparts[3])
#    resol=float(lineparts[4])
#    insertstr="Insert into structures values('%s','%s','%f','%f','%f')" % \
#    (pdb_id,chain_id,resol,rwork,rfree)
#    cursor.execute(insertstr)
#    db.commit()
#infile.close()
count=0
infile=open("list_of_maps_new_temp.txt","r")
for line in infile:
    count+=1
    lineparts=line.split("_")
    pdb_id=lineparts[0]
    chain_id=lineparts[1]
    if count%1000==0:
        print pdb_id+","+chain_id
    resname=lineparts[2][0:3]
    resno=int(lineparts[2][3:])
    molprobfile=open("ramalyze_outputs/"+pdb_id+".rama","r")
    x=molprobfile.readline()
    for line2 in molprobfile:
        line2parts=line2.split(":")
        resdetails=line2parts[0]
        mol_chain=resdetails[1]
        mol_resno=resdetails[2:7]
        if mol_resno[-1]!=" ":
            continue
        mol_resno=int(resdetails[2:7].strip())
        if mol_chain==chain_id and mol_resno==resno:
            phi=int(float(line2parts[2]))
            psi=int(float(line2parts[3]))
            break
    molprobfile.close()
    rotate_file=open("peptide_rotation_proteins/"+line[0:-10]+".csv","r")
    x=rotate_file.readline()
    for line1 in rotate_file:
        line1parts=line1.split(",")
        if phi==int(line1parts[0]) and psi==int(line1parts[1]):
            apd=int(line1parts[2][0:-1])
            break
    rotate_file.close()
    insertstr="Insert into residues values('%s','%s','%s','%d','%d','%d','%s')" % \
    (pdb_id,chain_id,resname,resno,phi,psi,apd)
    cursor.execute(insertstr)
    db.commit()
infile.close()
        