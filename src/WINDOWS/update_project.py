#! /usr/bin/env python
# Utility to update XML Visual Studio projects
# @author Richard Berger <richard.berger@jku.at>
import os
import re
import sys
import glob
from xml.dom.minidom import parseString
import xml.etree.ElementTree as ET

if len(sys.argv) < 2:
  print("USAGE python3 update_project.py [FILENAME.vcxproj]")
  sys.exit(-1)

filename = sys.argv[1]

if not os.path.exists(filename):
  print("File ", filename, " not found!")
  sys.exit(-1)

fileName, fileExtension = os.path.splitext(filename)

if fileExtension != '.vcxproj':
  print("*.vcxproj extension required! Only XML-based project files supported!")
  sys.exit(-1)

ET.register_namespace('', 'http://schemas.microsoft.com/developer/msbuild/2003')
tree = ET.parse(filename)
root = tree.getroot()

header_files = glob.glob(os.path.join('..', '*.h'))
impl_files   = glob.glob(os.path.join('..', '*.cpp'))

for child in root:
  if len(child) > 0:
    if 'ItemGroup' in child.tag and 'ClInclude' in child[0].tag:
      # update header files ItemGroup
      child.clear()

      # add erfc implementation (required by VS < 2013)
      ET.SubElement(child, 'ClInclude', {'Include' : "extra\erf.h"})

      for hfile in header_files:
        relpath = '..\\' + os.path.basename(hfile)
        ET.SubElement(child, 'ClInclude', {'Include' : relpath})

    elif 'ItemGroup' in child.tag and 'ClCompile' in child[0].tag:
      # update implementation files ItemGroup
      child.clear()

      # add erfc implementation (required by VS < 2013)
      ET.SubElement(child, 'ClCompile', {'Include' : "extra\erf_namd.c"})

      for ifile in impl_files:
        relpath = '..\\' + os.path.basename(ifile)
        ET.SubElement(child, 'ClCompile', {'Include' : relpath})

txt = ET.tostring(root)
txt = parseString(txt).toprettyxml()
txt = re.sub('\s*\n(\s*\n)*', '\n', txt)

with open(filename, 'w') as f:
  f.write(txt)
