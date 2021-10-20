#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Create a .c file which initializes all of the exters
in a .h file


Command line usage (if any):

    usage: init_extern.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

211019 ksl Coding begun

'''

import sys
import os
from astropy.io import ascii
import numpy as np
import shutil


python_header='''
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "atomic.h"

'''



def doit(filename='atomic.h',outputfile=''):
    '''
    Do something magnificent

    Description:

    Notes:

    We try to capture comments 

    We also leave one blank line where there
    were multiple blank lines in the .h file

    History:


    '''

    try:
        x=open(filename)
        lines=x.readlines()
    except:
        print('Could not find ',filename)
        return

    if outputfile=='':
        outputfile=filename.replace('.h','_extern_init.c')

    if os.path.exists(outputfile):
        xout=outputfile.replace('.c','.old.c')
        shutil.copyfile(outputfile,xout)


    g=open(outputfile,'w')

    if filename=='python.h':
        g.write(python_header)
    g.write('#include "%s" \n\n '  % filename)

    comment=False
    empty=False 
    for one in lines:
        xone=one.strip()
        words=xone.split()
        if len(words)==0:
            if empty==False:
                g.write('\n')
            empty=True
        elif comment==True:
            g.write(one)
            empty=False
            if one.count('*/'):
                comment=False
        elif words[0]=='extern':
            xone=one.replace('extern ','')
            g.write(xone)
            empty=False
            if xone.count('/*') and not xone.count('*/'):
                comment=True
            else:
                comment=False



            

    

    return


    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print ('usage: init_extern.py filename')
