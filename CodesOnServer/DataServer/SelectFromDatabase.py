#!/usr/bin/env python
""" This code is an example to select required data from TIRSPEC's SQL database of all the observation logs"""

DatabaseFilename = '~/TIRSPECDataLog.db'

import sqlite3

with sqlite3.connect(DatabaseFilename) as con:
    c = con.cursor()
    c.execute(" SELECT Directory FROM DirectoryTable")   # Look in all direcotries in DB
    Directories = [row[0] for row in c]
    for D_Dir in Directories:
        print('***'+D_Dir+'***')
        c.execute("SELECT * FROM {0:s} WHERE lower= ?".format(D_Dir),('G',)) # Look for Spectral data 
        for row in c:
            print(' '.join([str(r) for r in row]))
