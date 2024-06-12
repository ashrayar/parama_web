#!/usr/bin/python
import cgi, cgitb
cgitb.enable()
import os
import sys
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders

form = cgi.FieldStorage() 
toaddr=form.getvalue('toaddr')
pdb_id=form.getvalue('pdb_id')
chain_id=form.getvalue('chain_id')
db_presence=form.getvalue('db_presence')
print "Content-type:text/html\r\n\r\n"
print
print '<html>'
print '<head>'
print '<title>Bond-parameter specific Ramachandran map</title>'
print '<link rel="stylesheet" type = "text/css" href = "parama.css" media = " all" />'
print '</head>'
print '<body>'
print '<div class="header">'
print '<img src="hide_square.png" align="left" style="width:120px;height:120px;"/>'
print '<img src="IISc_logo.png" align="right" style="width:120px;height:120px;"/>'
print '	<h2 style="color:#3399ff;">BORAM</h2>'
print '<p style="font-size:25px;font-family: Verdana,Geneva, sans-serif; text-align:center;color:white;"><b style="color:#3399ff;">Bo</b>nd-parameter specific <b style="color:#3399ff;">Ra</b>machandran <b style="color:#3399ff;">M</b>aps </p>'
print '</div>'
print '<div class="topnav">'
print '	<a href="index.html">Home</a>'
print ' <a href="about.html">About BORAM</a>'
print '	<a href="help.html">Help</a>'
print '	<a href="publication.html">Publication</a>'
print '	<a href="contact.html" style="float:right">Contact Us</a>'
print '</div>'
print '<div class="row">'
#print '<h4> DB presence: '+db_presence+' </h4>'
print '<h4> You will receive a mail with the link to detailed results of '+pdb_id+' to '+toaddr+'. Thank You </h4>'
print '</div>'
print '</body>'
print '</html>'
sys.stdout.flush()
os.system("at now <<< 'python protein_pep_unit_rotation.py "+pdb_id+" "+chain_id+" "+db_presence+" "+toaddr+"'")
