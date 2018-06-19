#!/usr/bin/python3
"""
read_calculation.py by Aksyonov Dmitry, Skoltech, Moscow


"""
import sys
sys.path.extend(['/home/aksenov/www/siman'])
import cgitb, cgi 
from analysis import calc_redox
from classes import CalculationVasp

db_path = '/home/aksenov/www/CES/_aksenov/'

cgitb.enable()
print('Content-type: text/html\n\n')
form = cgi.FieldStorage() 

# Get data from fields
string = form.getvalue('param')
string = string.replace('"', '')
files = string.split()


cl1 = CalculationVasp().deserialize(db_path+files[0])
cl2 = CalculationVasp().deserialize(db_path+files[1])

# output = cl1.energy_sigma0
output = calc_redox(cl1, cl2)['redox_pot']

print(output, 'V')


# sys.stdout.write(str(output) )
# sys.stdout.flush()
sys.exit(0)



