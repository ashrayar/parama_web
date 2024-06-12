#!/usr/bin/python
import MySQLdb
import datetime
def connect_to_db():
	host="localhost"
	username="root"
	db_name="boram"
	password="mysql123"
	db = MySQLdb.connect(host,username,password,db_name)
	cursor=db.cursor()
	return cursor

	
def check_present_in_db(pdb_id,chain_id="not_entered"):
	db_cursor=connect_to_db()
	if chain_id=="not_entered":
		sql="select chain_id from residues where pdb_id='"+pdb_id+"';"
		db_cursor.execute(sql)
		results=db_cursor.fetchall()
		if len(results)==0:
			return "False"
		else:
			return results[0][0]
	else:
		sql="select resnumber from residues where pdb_id='"+pdb_id+"' and chain_id='"+chain_id+"';"
		db_cursor.execute(sql)
		results=db_cursor.fetchall()
		if len(results)==0:
			return "False"
		else:
			return "True"
	
def get_results_from_db(pdb_id,chain_id,residues=[]):
	logfile=open("rotation_log.log","a")
	logfile.write("Inside get_results_from_db\n")
	logfile.close()
	db_cursor=connect_to_db()
	if len(residues)==0:
		sql="select * from residues where pdb_id='"+pdb_id+"' and chain_id='"+chain_id+"' order by resnumber;"
	else:
		reslist="("
		for item in residues:
			reslist+=(str(item)+",")
		reslist=reslist[0:-1]
		reslist+=(")")
		sql="select * from residues where pdb_id='"+pdb_id+"' and chain_id='"+chain_id+"' and resnumber in "+reslist+" order by resnumber;"
	db_cursor.execute(sql)
	results=db_cursor.fetchall()
	if len(results)==0:
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
		htmlfile.write('<h2 style="color:#3399ff;"></h2>\n')
		htmlfile.write('<p style="font-size:25px;font-family: Verdana,Geneva, sans-serif; text-align:center;color:white;"><b style="color:#3399ff;"></p>\n')
		htmlfile.write('</div>\n')
		htmlfile.write('<h3> All of the entered residues either were missing entirely or missing some backbone atoms (including hydrogens).</h3>')
		htmlfile.close()
	else:
		cur_date=datetime.datetime.now()
		uniq_id=str(cur_date.year)+str(cur_date.month)+str(cur_date.day)+"_"+str(cur_date.hour)+str(cur_date.minute)+str(cur_date.second)
		htmlfile=open('results/'+pdb_id+'_'+chain_id+'_'+uniq_id+'.html','w')
		#htmlfile.write('<html>\n')
		#htmlfile.write('<head>\n')
		#htmlfile.write('<title>Bond-parameter specific Ramchandran Maps</title>\n')
		#htmlfile.write('<link rel="stylesheet" type = "text/css" href = "parama.css" media = " all" />\n')
		#htmlfile.write('<style>\n')
		#htmlfile.write('td { border: 1px solid #ddd; padding:8px; background-color: #f2f2f2; text-align:center;}\n')
		#htmlfile.write('th {background-color: #3399ff;border: 1px solid black;padding:8px;text-align:center;}\n')
		#htmlfile.write('tr:nth-child(even) {background-color: #f2f2f2;}\n')
		#htmlfile.write('</style>\n')
		#htmlfile.write('</head>\n')
		#htmlfile.write('<body>\n')
		#htmlfile.write('<div class="header">\n')
		#htmlfile.write('<h2 style="color:#3399ff;"></h2>\n')
		#htmlfile.write('<p style="font-size:25px;font-family: Verdana,Geneva, sans-serif; text-align:center;color:white;"><b style="color:#3399ff;"></p>\n')
		htmlfile.write('</div>\n')
		htmlfile.write('<style>\n')
		htmlfile.write('td { border: 1px solid #ddd; padding:8px; background-color: #f2f2f2; text-align:center;}\n')
		htmlfile.write('th {background-color: #3399ff;border: 1px solid black;padding:8px;text-align:center;}\n')
		htmlfile.write('tr:nth-child(even) {background-color: #f2f2f2;}\n')
		htmlfile.write('</style>\n')
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
		for row in results:
			htmlfile.write('<tr>\n')
			htmlfile.write('<td>'+str(row[3])+'</td>\n')
			htmlfile.write('<td>'+row[2]+'</td>\n')
			map_image_file="http://pauling.mbu.iisc.ac.in/bond_param_rmap/peptide_rmaps_proteins/"+row[0]+"_"+row[1]+"_"+row[2]+str(row[3])+"_rmap.png"
			htmlfile.write('<td><a href="'+map_image_file+'">Image file</a></td>\n')
			rotation_file="http://pauling.mbu.iisc.ac.in/bond_param_rmap/peptide_rotation_proteins/"+row[0]+"_"+row[1]+"_"+row[2]+str(row[3])+".csv"
			htmlfile.write('<td><a href="'+rotation_file+'">Data File</a></td>\n')
			htmlfile.write('<td>'+str(row[4])+'</td>\n')
			htmlfile.write('<td>'+str(row[5])+'</td>\n')
			if row[6]=="1":
				htmlfile.write('<td>This residue\'s &straightphi;-&psi; is fully allowed according to its Ramachandran Map</td>\n')
			elif row[6]=="2":
				htmlfile.write('<td>This residue\'s &straightphi;-&psi; is partially allowed according to its Ramachandran Map</td>\n')
			else:
				htmlfile.write('<td>This residue\'s &straightphi;-&psi; is disallowed according to its Ramachandran Map</td>\n')
			htmlfile.write('</tr>\n')
		htmlfile.write('</table>\n')
		htmlfile.write('</body>\n')
		htmlfile.write('</html>\n')
		htmlfile.close()
	return uniq_id

#print check_present_in_db("6t54")	


			
			
		
