#!/usr/bin/python

import math
import numpy as np
from plot_generated_rmaps import plot_rmap
def FindDihedralAngle(A,B,C,D):
    BCi=C[0]-B[0]
    BCj=C[1]-B[1]
    BCk=C[2]-B[2]

    BAi=A[0]-B[0]
    BAj=A[1]-B[1]
    BAk=A[2]-B[2]

    CDi=D[0]-C[0]
    CDj=D[1]-C[1]
    CDk=D[2]-C[2]

    Q1i=(BCj*BAk)-(BCk*BAj)
    Q1j=(BCk*BAi)-(BCi*BAk)
    Q1k=(BCi*BAj)-(BCj*BAi)

    Q2i=(BCj*CDk)-(BCk*CDj)
    Q2j=(BCk*CDi)-(BCi*CDk)
    Q2k=(BCi*CDj)-(BCj*CDi)
    magQ1=math.sqrt((Q1i*Q1i)+(Q1j*Q1j)+(Q1k*Q1k))
    Q1i=Q1i/magQ1
    Q1j=Q1j/magQ1
    Q1k=Q1k/magQ1

    magQ2=math.sqrt((Q2i*Q2i)+(Q2j*Q2j)+(Q2k*Q2k))
    Q2i=Q2i/magQ2
    Q2j=Q2j/magQ2
    Q2k=Q2k/magQ2

    Q1dotQ2=(Q1i*Q2i)+(Q1j*Q2j)+(Q1k*Q2k)
    chi=math.acos(Q1dotQ2)
    chinew=math.degrees(chi)
    
    Q1=np.array([Q1i,Q1j,Q1k])
    Q2=np.array([Q2i,Q2j,Q2k])
    Q1crossQ2=np.cross(Q1,Q2)
    magBC=math.sqrt((BCi*BCi)+(BCj*BCj)+(BCk*BCk))
    unitBCi=BCi/magBC
    unitBCj=BCj/magBC
    unitBCk=BCk/magBC
    unitBC=np.array([unitBCi,unitBCj,unitBCk])
    anglesign=np.dot(Q1crossQ2,unitBC)
    if anglesign<0:
        chinew=chinew*-1
        
    return int(round(chinew))

def resultantpoint(l,m,n,a,b,c,x,y,z,angle):
    angle=math.radians(angle)
    term1=(a*(m*m+n*n)-l*(b*m+c*n-l*x-m*y-n*z))*(1-math.cos(angle))
    term2=x*math.cos(angle)
    term3=(b*n-c*m-n*y+m*z)*math.sin(angle)
    x_new=term1+term2+term3
    term1=(b*(l*l+n*n)-m*(a*l+c*n-l*x-m*y-n*z))*(1-math.cos(angle))
    term2=y*math.cos(angle)
    term3=(c*l-a*n+n*x-l*z)*math.sin(angle)
    y_new=term1+term2+term3
    term1=(c*(l*l+m*m)-n*(a*l+b*m-l*x-m*y-n*z))*(1-math.cos(angle))
    term2=z*math.cos(angle)
    term3=(a*m-b*l-m*x+l*y)*math.sin(angle)
    z_new=term1+term2+term3
    newcoords=[0 for z in range(3)]
    newcoords[0]=x_new
    newcoords[1]=y_new
    newcoords[2]=z_new
    return newcoords

def FindDistance(a1,a2):
    return round(math.sqrt(((a2[0]-a1[0])*(a2[0]-a1[0]))+((a2[1]-a1[1])*(a2[1]-a1[1]))+((a2[2]-a1[2])*(a2[2]-a1[2]))),3)

def CheckShortContacts(a1,a2,d1,d2):
    dist=FindDistance(a1,a2)
    if dist<=d1 and dist>=d2:
        return "partially allowed: "+str(dist)
    elif dist>d1:
        return "fully allowed: "+str(dist)
    elif dist<d2:
        return "disallowed: "+str(dist)
    

def BringToZero(tempatoms,prepro,tempphi,temppsi):
    #bring phi to 0
    dircos=FindDirCosines(tempatoms['N2'], tempatoms['CA2'])
    tempatoms['CB']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CB'][0],tempatoms['CB'][1],tempatoms['CB'][2],-tempphi)
    tempatoms['CA2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],-tempphi)
    tempatoms['HA']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['HA'][0],tempatoms['HA'][1],tempatoms['HA'][2],-tempphi)
    tempatoms['C2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['C2'][0],tempatoms['C2'][1],tempatoms['C2'][2],-tempphi)
    tempatoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['O2'][0],tempatoms['O2'][1],tempatoms['O2'][2],-tempphi)
    tempatoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['N3'][0],tempatoms['N3'][1],tempatoms['N3'][2],-tempphi)
    if prepro=="normal":
        tempatoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['H3'][0],tempatoms['H3'][1],tempatoms['H3'][2],-tempphi)
    else:
        tempatoms['CD3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CD3'][0],tempatoms['CD3'][1],tempatoms['CD3'][2],-tempphi)
    tempatoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CA3'][0],tempatoms['CA3'][1],tempatoms['CA3'][2],-tempphi)
    #bring psi to 0
    dircos=FindDirCosines(tempatoms['CA2'],tempatoms['C2'])
    tempatoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['O2'][0],tempatoms['O2'][1],tempatoms['O2'][2],-temppsi)
    tempatoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['N3'][0],tempatoms['N3'][1],tempatoms['N3'][2],-temppsi)
    if prepro=="normal":
        tempatoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['H3'][0],tempatoms['H3'][1],tempatoms['H3'][2],-temppsi)
    else:
        tempatoms['CD3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['CD3'][0],tempatoms['CD3'][1],tempatoms['CD3'][2],-temppsi)
    tempatoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['CA3'][0],tempatoms['CA3'][1],tempatoms['CA3'][2],-temppsi)

def add_missing_phipsi(phipsi):
    keys=phipsi.keys()
    if "-178,-180" not in keys and "-178,180" in keys:
        for i in range(-179,180):
            phipsi[str(i)+","+"-180"]=phipsi[str(i)+","+"180"]
    elif "-178,180" not in keys and "-178,-180" in keys:
        for i in range(-179,180):
            phipsi[str(i)+","+"180"]=phipsi[str(i)+","+"-180"]
    if "-180,-178" not in keys and "180,-178" in keys:
        for j in range(-179,180):
            phipsi["-180,"+str(j)]=phipsi["180,"+str(j)]
    elif "180,-178" not in keys and "-180,-178" in keys:
        for j in range(-179,180):
            phipsi["180,"+str(j)]=phipsi["-180,"+str(j)]

    if "-180,180" in keys:
        phipsi["180,180"]=phipsi["-180,180"]
        phipsi["180,-180"]=phipsi["-180,180"]
        phipsi["-180,-180"]=phipsi["-180,180"]
    elif "-180,-180" in keys:
        phipsi["180,180"]=phipsi["-180,-180"]
        phipsi["180,-180"]=phipsi["-180,-180"]
        phipsi["-180,180"]=phipsi["-180,-180"]
    elif "180,-180" in keys:
        phipsi["180,180"]=phipsi["180,-180"]
        phipsi["-180,180"]=phipsi["180,-180"]
        phipsi["-180,-180"]=phipsi["180,-180"]
    elif "180,180" in keys:
        phipsi["-180,180"]=phipsi["180,180"]
        phipsi["180,-180"]=phipsi["180,180"]
        phipsi["-180,-180"]=phipsi["180,180"]


    
def Ala_MapsGenerator(atoms,pdb_id,chain_id,resname,resnumber,prepro):
    digimapfile=open("results/"+pdb_id+"_"+chain_id+"_"+resname+str(resnumber)+".csv","w")
    digimapfile.write("phi,psi,a/d/p\n")
    actphi=FindDihedralAngle(atoms['C1'],atoms['N2'],atoms['CA2'],atoms['C2'])
    actpsi=FindDihedralAngle(atoms['N2'],atoms['CA2'],atoms['C2'],atoms['N3'])
    BringToZero(atoms,prepro,actphi,actpsi)
    phipsi=dict()
    for i in range(-180,180):
        for j in range(-180,180):
            finalres=''
            #bring psi to 0
            dircos=FindDirCosines(atoms['CA2'],atoms['C2'])
            atoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['O2'][0],atoms['O2'][1],atoms['O2'][2],1)
            atoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['N3'][0],atoms['N3'][1],atoms['N3'][2],1)
            if prepro=="normal":
                atoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['H3'][0],atoms['H3'][1],atoms['H3'][2],1)
            else:
                atoms['CD3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['CD3'][0],atoms['CD3'][1],atoms['CD3'][2],1)
            atoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['CA3'][0],atoms['CA3'][1],atoms['CA3'][2],1)
            phi=FindDihedralAngle(atoms['C1'],atoms['N2'],atoms['CA2'],atoms['C2'])
            psi=FindDihedralAngle(atoms['N2'],atoms['CA2'],atoms['C2'],atoms['N3'])
            #sampleout.writelines("\nphi="+str(phi)+"  psi="+str(psi)+"\n")
            if prepro=="normal":            
                res=CheckShortContacts(atoms['CA1'],atoms['H3'],2.40,2.20)
            else:
                res=CheckShortContacts(atoms['CA1'],atoms['CD3'],3.00,2.90)
            finalres=finalres+res[0]
            #if res[0]=='p' or res[0]=='d':
                #sampleout.write('CA1...H3 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['HA'],2.40,2.20)
            finalres=finalres+res[0]
            #if res[0]=='p' or res[0]=='d':
                #sampleout.write('C1...HA : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['C2'],3.00,2.90)
            finalres=finalres+res[0]
            #if res[0]=='p' or res[0]=='d':
                #sampleout.write('C1...C2 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['O2'],2.80,2.70)
            finalres=finalres+res[0]
            #if res[0]=='p' or res[0]=='d':
                #sampleout.write('C1...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['N3'],2.90,2.80)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('C1...N3 : '+res+'\n')
            if prepro=="normal":  
                res=CheckShortContacts(atoms['C1'],atoms['H3'],2.40,2.20)
            else:
                res=CheckShortContacts(atoms['C1'],atoms['CD3'],3.00,2.90)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('C1...H3 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['CB'],3.20,3.00)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('C1...CB : '+res+'\n')
            
            res=CheckShortContacts(atoms['O1'],atoms['CB'],2.80,2.70)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...CB : '+res+'\n')
            res=CheckShortContacts(atoms['O1'],atoms['HA'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...HA : '+res+'\n')
            res=CheckShortContacts(atoms['O1'],atoms['C2'],2.80,2.70)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...C2 : '+res+'\n')    
            res=CheckShortContacts(atoms['O1'],atoms['O2'],2.70,2.60)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...O2 : '+res+'\n')  
            res=CheckShortContacts(atoms['O1'],atoms['N3'],2.70,2.60)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...N3 : '+res+'\n') 
            if prepro=="normal":  
                res=CheckShortContacts(atoms['O1'],atoms['H3'],2.40,2.20)
            else:
                res=CheckShortContacts(atoms['O1'],atoms['CD3'],2.80,2.70)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...H3 : '+res+'\n') 
            res=CheckShortContacts(atoms['O1'],atoms['CA3'],2.80,2.70)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...CA3 : '+res+'\n')  

            res=CheckShortContacts(atoms['N2'],atoms['O2'],2.70,2.60)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('N2...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['N2'],atoms['N3'],2.70,2.60)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('N2...N3 : '+res+'\n')  
            if prepro=="normal":  
                res=CheckShortContacts(atoms['N2'],atoms['H3'],2.40,2.20)
            else:
                res=CheckShortContacts(atoms['N2'],atoms['CD3'],2.90,2.80)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('N2...H3 : '+res+'\n')
            
            res=CheckShortContacts(atoms['H2'],atoms['HA'],2.00,1.90)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...HA : '+res+'\n')   
            res=CheckShortContacts(atoms['H2'],atoms['C2'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...C2 : '+res+'\n')
            res=CheckShortContacts(atoms['H2'],atoms['O2'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['H2'],atoms['N3'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...N3 : '+res+'\n')
            if prepro=="normal":  
                res=CheckShortContacts(atoms['H2'],atoms['H3'],2.00,1.90)
            else:
                res=CheckShortContacts(atoms['H2'],atoms['CD3'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...H3 : '+res+'\n')
            res=CheckShortContacts(atoms['H2'],atoms['CB'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...CB : '+res+'\n')
            
            res=CheckShortContacts(atoms['CB'],atoms['O2'],2.80,2.70)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('CB...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['CB'],atoms['N3'],2.90,2.80)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('CB...N3 : '+res+'\n')
            if prepro=="normal":  
                res=CheckShortContacts(atoms['CB'],atoms['H3'],2.40,2.20)
            else:
                res=CheckShortContacts(atoms['CB'],atoms['CD3'],3.20,3.00)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('CB...H3 : '+res+'\n')
            
            res=CheckShortContacts(atoms['HA'],atoms['O2'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('HA...O2 : '+res+'\n')
            
            res=CheckShortContacts(atoms['HA'],atoms['N3'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('HA...N3 : '+res+'\n')
            if prepro=="normal":  
                res=CheckShortContacts(atoms['HA'],atoms['H3'],2.00,1.90)
            else:
                res=CheckShortContacts(atoms['HA'],atoms['CD3'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('HA...H3 : '+res+'\n')
            
            if 'd' not in finalres and 'p' not in finalres:
#                sampleout.write("This conformation is fully allowed\n--------------------------------------------------------\n")
#                digimapfile.write(str(phi)+","+str(psi)+",1\n")
                phipsi[str(phi)+","+str(psi)]=1
            elif 'd' in finalres:
#                sampleout.write("This conformation is disallowed\n--------------------------------------------------------\n")
#                digimapfile.write(str(phi)+","+str(psi)+",0\n")
                phipsi[str(phi)+","+str(psi)]=0
            else:
#                sampleout.write("This conformation is partially allowed\n--------------------------------------------------------\n")
#                digimapfile.write(str(phi)+","+str(psi)+",2\n") 
                phipsi[str(phi)+","+str(psi)]=2
        
        
        dircos=FindDirCosines(atoms['N2'], atoms['CA2'])
        atoms['CB']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CB'][0],atoms['CB'][1],atoms['CB'][2],1)
        atoms['CA2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],1)
        atoms['HA']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['HA'][0],atoms['HA'][1],atoms['HA'][2],1)
        atoms['C2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['C2'][0],atoms['C2'][1],atoms['C2'][2],1)
        atoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['O2'][0],atoms['O2'][1],atoms['O2'][2],1)
        atoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['N3'][0],atoms['N3'][1],atoms['N3'][2],1)
        if prepro=="normal":
            atoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['H3'][0],atoms['H3'][1],atoms['H3'][2],1)
        else:
            atoms['CD3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CD3'][0],atoms['CD3'][1],atoms['CD3'][2],1)
        atoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CA3'][0],atoms['CA3'][1],atoms['CA3'][2],1)   
            
            
                    
    #sampleout.close()
    #digimapfile.write(str(actphi)+","+str(actpsi)+",3")
    add_missing_phipsi(phipsi)
    fa_phi=[]
    fa_psi=[]
    pa_phi=[]
    pa_psi=[]
    for i in range(-180,181):
        for j in range(-180,181):
            digimapfile.write(str(i)+","+str(j)+","+str(phipsi[str(i)+","+str(j)])+"\n")
            if phipsi[str(i)+","+str(j)]==1:
                fa_phi.append(i)
                fa_psi.append(j)
            elif phipsi[str(i)+","+str(j)]==2:
                pa_phi.append(i)
                pa_psi.append(j)
    digimapfile.close()
    plot_rmap(fa_phi,fa_psi,pa_phi,pa_psi,actphi,actpsi,pdb_id,chain_id,resnumber,resname)
    return [actphi,actpsi,phipsi[str(actphi)+","+str(actpsi)]]
    
def FindDirCosines(A,B):
    delx=B[0]-A[0]
    dely=B[1]-A[1]
    delz=B[2]-A[2]
    denom=math.sqrt(delx*delx+dely*dely+delz*delz)
    l=delx/denom
    m=dely/denom
    n=delz/denom
    a=[l,m,n]
    return a

    
        
            
            
        
        
