# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:29:30 2026

Helper module for the work with amdp from a Python IDE.
It calls amdp with the required argument(s)

@author: Laszlo Merenyi
"""

import sys
import runpy

# The last uncommented project directory will be used

project_path = r"d:\DIM-data\projetc_1"
project_path = r"d:\DIM-data\projetc_2"


sys.argv = ['amdp.py', project_path]

runpy.run_path('amdp.py', run_name='__main__')


     
     