# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:29:30 2026

Helper module for the work with amdp from a Python IDE.
It calls amdp with the required argument(s)

@author: Laszlo Merenyi
"""

import sys
import runpy

project_path = r"d:\DIM-data\projetc_1"
project_path = r"d:\DIM-data\projetc_2"

# --------- Remove these from public release --------------
# project_path = r"A:\ADS\BFE\ADS-5 2026-06"
# project_path = r"A:\ADS\MNY\ADS-7 2026\2026-06-13"
# project_path = r"A:\ADS\MNY\ADS-8 2026\2026-06-19"
project_path = r"A:\ADS\MNY\ADS-9 2026\2026-07-30"

# project_path = r"P:\dev\datacq\ADS\amdp\configs\NAQ+BFE\BFE project"
# project_path = r"P:\dev\datacq\ADS\amdp\configs\NAQ+BFE\NAQ project"

# project_path = r"P:\dev\datacq\ADS\amdp\configs\MNY\A1-act"
# project_path = r"P:\dev\datacq\ADS\amdp\configs\MNY\A2-act"
# project_path = r"A:\ADS\SOD\ADS3-2025-automn"
# project_path = r"P:\dev\datacq\ADS\amdp\configs\MNY-absmeas"


# project_path = r"A:\ADS\_ACT\MNY-A1-parameters" 
project_path = r"A:\ADS\_ACT\MNY-A2-parameters" 
project_path = r"D:\support\Columbia-timestamp-problem"

# --------- End of remove --------------

sys.argv = ['amdp.py', project_path]

runpy.run_path('amdp.py', run_name='__main__')


     
     