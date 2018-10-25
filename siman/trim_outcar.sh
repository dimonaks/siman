#!/bin/bash
nbands=`grep "NBANDS=" filename | awk '{{print $NF - 1}}'`; 
sed -i -e "/band No./,+${{nbands}}d" filename