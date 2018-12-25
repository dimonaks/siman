#!/usr/bin/python

import os, sys
import xml.etree.ElementTree as ET


root = ET.fromstring(sys.stdin.read())

# for child in root:
#   print( child.tag, child.attrib)

for job in root.findall('Job'):
    # name = job.get('Job_Name')
    name = job.find('Job_Name').text
    path = job.find('Output_Path').text
    path = os.path.dirname(path)
    path = path.replace('arcuda:', 'cd ')
    print(path)
