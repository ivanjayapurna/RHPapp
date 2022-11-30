# Author: Spencer Jenkins
# Contains file directories and other global variables.

import os

XUGROUP_IGOR_MACRO_SIZE = 39

# These directories apply in development mode.
OUTPUTS_DIR = "app/static/outputs/"
PLOTS_DIR = 'app/static/plots/'
# This directory should be changed to the offline data directory for use in development mode.
OFFLINE_DATA_DIR = "F:\\summer-2021-polymers\\data\\"

OCF = True  # CHANGE TO FALSE IF NOT ON OCF
OCF_USERNAME = 'spencerrjenkins'  # SET TO OCF USERNAME

# SQL Database parameters.
OCF_USERNAME = 'xugroup'  # SET TO OCF USERNAME
DB_HOST = 'mysql'
DB_USERNAME = OCF_USERNAME
DB_NAME = OCF_USERNAME
DB_PASSWORD = 'b41HbOoxcNIfswr87STGtN22'

if OCF:
    OUTPUTS_DIR = os.path.expanduser('~/public_html') + '/static/outputs/'
    PLOTS_DIR = os.path.expanduser('~/public_html') + '/static/plots/'

HLB_CUTOFF = 9

