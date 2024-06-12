#!/usr/bin/python
import cgi, cgitb
cgitb.enable()
import os
import urllib
from test_phipsi_pdb import check_chain_id
from test_phipsi_pdb import addchainid
from DbOperations import check_present_in_db
from DbOperations import get_results_from_db
form = cgi.FieldStorage() 
print "Content-type:text/html\r\n\r\n"
print
print '<html>'
print '<head>'
print '<title>Bond-parameter specific Ramachandran maps</title>'
print '<link rel="stylesheet" type = "text/css" href = "parama.css" media = " all" />'
print '</head>'
print '<body>'
print '<div class="header">'
print '<img src="hide_square.png" align="left" style="width:120px;height:120px;"/>'
print '<img src="IISc_logo.png" align="right" style="width:120px;height:120px;"/>'
print '<h2 style="color:#3399ff;">BORAM</h2>'
print '<p style="font-size:25px;font-family: Verdana,Geneva, sans-serif; text-align:center;color:white;"><b style="color:#3399ff;">Bo</b>nd-parameter specific <b style="color:#3399ff;">Ra</b>machandran <b style="color:#3399ff;">M</b>aps </p>'
print '</div>'
flag=0
ready_to_submit=0
UPLOAD_DIR = 'pdb_files/'
# Get data from fields
fileitem=form['pdb_file']
#check if chain ID has been entered
if form.getvalue('chain_id'):
	chain_id=form.getvalue('chain_id').strip().upper()
else:
	chain_id="not_entered"
if form.getvalue('pdb_id'): #pdb id entered, try to download
	pdb_id=form.getvalue('pdb_id')
	if chain_id=="not_entered":
		db_presence=check_present_in_db(pdb_id)
		if db_presence=="False":
			present_in_db=False
		else:
			present_in_db=True
			chain_id=db_presence
	else:
		db_presence=check_present_in_db(pdb_id,chain_id)
		if db_presence=="False":
			present_in_db=False
		else:
			present_in_db=True
	if not present_in_db:
		url="https://files.rcsb.org/download/"+pdb_id+".pdb"
		uh=urllib.urlopen(url)
		file_contents=uh.read()
		if "<!DOCTYPE HTML" in file_contents: #pdb id does not exist, returns an html error file
			print '<h3> Invalid PDB ID </h3>'
			print '</div>'
			print '</body>'
			print '</html>'
			flag=1
		else: #valid pdb id downloaded, save locally
			pdbfile=open('pdb_files/'+pdb_id+'.pdb','w')
			pdbfile.write(file_contents)
			pdbfile.close()
			pdb_id=pdb_id+".pdb"
			pdbfile=open('pdb_files/'+pdb_id,'r') #open pdb file and keep
			pdb_id=pdb_id[0:-4]
			#print '<h3> PDB file has been downloaded from Protein Data Bank </h3>'
	else:
		ready_to_submit=1
elif fileitem.filename: #pdb file has been uploaded
	pdb_id=fileitem.filename
	pdb_id=pdb_id[0:-4]
	if chain_id=="not_entered":
		db_presence=check_present_in_db(pdb_id)
		if db_presence=="False":
			present_in_db=False
		else:
			present_in_db=True
			chain_id=db_presence
	else:
		db_presence=check_present_in_db(pdb_id,chain_id)
		if db_presence=="False":
			present_in_db=False
		else:
			present_in_db=True
	if not present_in_db:
		uploaded_file_path = os.path.join(UPLOAD_DIR, os.path.basename(fileitem.filename))
		with file(uploaded_file_path, 'wb') as fout: #save the uploaded file locally
			while True:
				chunk = fileitem.file.read(100000)
				if not chunk:
					break
				fout.write (chunk)
		fout.close()
		#print '<h3> PDB file has been uploaded to server</h3>'
	else:
		ready_to_submit=1
if flag!=1: #flag = 1 when invalid pdb id is entered
	if chain_id=='not_entered':
 		if not present_in_db:
			chain_id=addchainid(pdb_id)
	if not present_in_db:
		result=check_chain_id(pdb_id,chain_id) #checks if the entered chain id exists in the PDB file
		if result==0: #chain not found
			print '<h3> Chain not found in PDB file </h3>'
			print '</div>'
			print '</body>'
			print '</html>'
		else:
			ready_to_submit=1
	if ready_to_submit==1: #chain found, generate maps
		if present_in_db:
			if form.getvalue('residues'):
				residues=form.getvalue('residues')
				wanted_res_str=residues.split(",")
				wanted_residues=[int(x) for x in wanted_res_str]
				wanted_residues.sort()
				uniq_id=get_results_from_db(pdb_id,chain_id,wanted_residues)
			else:
				uniq_id=get_results_from_db(pdb_id,chain_id)
			htmlfile=open("results/"+pdb_id+"_"+chain_id+"_"+uniq_id+".html")
			print '<h3>The structure you have entered is available in our database</h3>'
			for line in htmlfile:
				print line
			htmlfile.close()
		else:
			print '<h3>The PDB file you have entered/uploaded is not available in our database. The results will be generated and sent to you via email</h3>'
			print '<script type="text/javascript">'
			print 'function validateEmail(){'
			print 'var re = /^(([^<>()\[\]\\.,;:\s@"]+(\.[^<>()\[\]\\.,;:\s@"]+)*)|(".+"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/;'
			print 'var email=document.getElementById("toaddr").value;'
			print 'if (re.test(email)){'
			print 'return true;'
			print '}'
			print 'alert("Please enter valid email address");'
			print 'return false;'
			print '}'
			print '</script>'
			if form.getvalue('residues'): # residue-wise maps need to be generated
				print '<form class="form1" method="post" action="ResidueMap_Interface.py" style="margin-left:50px;" onsubmit="return validateEmail();">'
				print '<input type="hidden" name="residues" value='+form.getvalue('residues')+'>'
				print '<input type="hidden" name="db_presence" value='+str(present_in_db)+'>'
			else: # generate maps for all residues in the given chain
				print '<form class="form1" method="post" action="MapGenerator_interface.py" style="margin-left:50px;" onsubmit="return validateEmail();">'
				print '<input type="hidden" name="db_presence" value='+str(present_in_db)+'>'
			if not present_in_db:	
				print 'Enter Email: <input type="text" name="toaddr" id="toaddr"><br><br>'
				print '<input type="hidden" name="pdb_id" value='+pdb_id+'>'
				print '<input type="hidden" name="chain_id" value='+chain_id+'>'
				print '<input class="button" type="submit" value="Submit">'
				print '</form>'	
	
print '</body>'
print '</html>'

