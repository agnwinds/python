#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Execute a series of short python runs to test whether a version 
of python is working.


Command line usage (if any):

    usage: regression.py [-np 3 -pf_dir test -out_dir foo] version 

    where 

        version         the executable of python
        -np 3           the number of processors with which to run (default 3)
        -pf_dir test    the directory containing all of the .pf files which will be run
                        The defaults is $PYTHON/examples/regress.  One does not need
                        to provide the full path name to the directory.  The routine doit
                        first searches the current workind directory for the directory and then
                        looks in $PYTHON/examples/
        -out_dir foo    The directory (below the current working directory) where the 
                        tests will run.  The defauld is constructed for the version
                        and the data

Description:  

    The basic process is as follows

    1.  Create a directory in which to work and intialize it
    2.  Copy all of the relevant .pf files to this directory
    3.  Switch to the working directory and run all the models, performing
        some basic checks 

Primary routines:

    doit:       Internal routine which runs python on all of the pf files of interest.  Use
                this if working in a python shell
    steer:      A routine to parse the command lineh
    check_one:  A routine which oversees checking of the runs

Notes:

    Regression here means to run a series of models. These routines do not compare the
    models to earlier runs

    It is possible to add special tests for regression, either to because models must be 
    run in a certain order or because special options must be used.  These should
    be added in a special section of doit, and should be self-contained.  

                                       
History:

170903 ksl Coding begun
180225 ksl nsh had added a special test for hydro.  This is useful in the sense that
           it allows one to deal with situations where several runs must be carreid out
           sequentially but but the mechanism that he chose makes
           it difficult to add new tests, because the standard and speciial tests are 
           not sufficently isolated from one another.


'''

import sys
from astropy.io import ascii
import numpy
from glob import glob
import time
import subprocess
import os
import shutil
import py_error
import balmer_decrement 
import regression_check




def sum_errors(root='1d_sn'):
    '''
    This sums the errors from various threads

    Note that this effectively duplicates a routine
    py_error which James wrote, but in Python 3


    '''
    diag_files=glob('diag_%s/%s_*.diag' % (root,root))

    errors=[]
    for one in diag_files:
        x=open(one)
        lines=x.readlines()
        lines=lines[len(lines)-100:]
        test=False
        for one_line in lines:
            if test and one_line.count(' -- '):
                errors.append(one_line)
            if one_line.count('Recurrences --  Description'):
                test=True

    # At this point we have the error but we need to parse them etc

    i=0
    error_names=[]
    error_counts=[]
    while i<len(errors):
        errors[i]=errors[i].strip()
        words=errors[i].split(' -- ')
        error_names.append(words[1])
        try:
            error_counts.append(int(words[0]))
        except ValueError:
            print('Error:sum_errors: %s in root %s' % (errors[i],root))
            error_counts.append(-999)
        i+=1

    # print(errors)

    # Now get the unique names

    unique_errors=set(error_names)

    # print(unique_errors)

    # finally add up the errors from the various threads

    records=[]
    for one in unique_errors:
        name=one
        counts=0
        i=0
        while i<len(error_names):
            if name==error_names[i]:
                counts+=error_counts[i]
            i+=1
        records.append([name,counts])
    return(records)


def check_one(xfile,root):
    '''
    Perform some checks on what happened in a run
    '''

    g=open(root+'.sig')
    lines=g.readlines()
    last_line=lines[len(lines)-1]
    words=last_line.split()
    if last_line.count('COMPLETE'):
        string='The model ran to completion in %s s' % (words[5])
    else:
        string='The models stopped before completing'

    print(string)
    xfile.write('%s\n' % string)

    errors=sum_errors(root)

    if len(errors)>0:
        xfile.write('\nThe reported errors in all threads was as follows: \n')
        for one in errors:
            xfile.write('%10d -- %s\n' % (one[1],one[0]))
    else:
        xfile.write('No diagnostic errors were reported\n')
    py_error.doit(root)
    return
	
	    
    


def doit(version='py',pf_dir='',out_dir='',np=3,outputfile='Summary.txt'):
    '''
    Test a specific version of python against a series of models

    Description:

    Notes:

        The routine looks for the input directory first as a subdirectory
        of the directory from which regression is run, and then in $PYTHON/examples/

    History:

        170903  ksl Bagan work
        170904  ksl Updated how routine looks for input directories
                    and attempted to get a more readable stderr output
                    Also, eliminted pf files with the extension .out.pf
                    because these are mostlikely duplicates

    '''

    date=time.strftime("%y%m%d", time.gmtime())
    

    if out_dir=='':
        out_dir='%s_%s' % (version,date)

    if os.path.exists(out_dir)==False:
        os.mkdir(out_dir)

    print(date)

    # Get the PYTHON environment variable

    PYTHON=os.environ['PYTHON']
    print(PYTHON)

    if pf_dir=='':
        pf_dir=PYTHON+'/examples/regress'
		
    if os.path.isdir(pf_dir):
        pf_files=glob(pf_dir+'/*pf')
        txt_files=glob(pf_dir+'/*.txt')
    elif os.path.isdir('%s/examples/%s' % (PYTHON,pf_dir)):
        pf_files=glob('%s/examples/%s/*pf' % (PYTHON,pf_dir))
        txt_files=glob('%s/examples/%s/*.txt' % (PYTHON,pf_dir))
    else:
        print('Error: The pf directory %s does not appear to exist' % pf_dir)
        return
		
    # screen out the .out.pf files because we may be running this on a directory
    # in which python has already been done

    select=[]
    for one in pf_files:
        if one.count('.out.pf')==0:
            select.append(one)
    pf_files=sorted(select)  # Assure that files are executed in alphabeticall order

    if len(pf_files)==0:
        print ('No input files found for %s search' % (pf_dir+'/*pf'))
        return
    else: 
        print ('Executing models for the following models:')
        for one in pf_files:
            print(one)



    commands=[]
    root_names=[]
    for one in pf_files:
        print('Copy',one)
        shutil.copy(one,out_dir)
        words=one.split('/')
        pf=(words[len(words)-1])
        root_name=pf.replace('.pf','')
        if np<=1:
            command='%s %s' % (version,pf)
        else:
            command='mpirun -np %d %s %s >%s.stdout.txt' % (np,version,pf,root_name)
        commands.append(command)
        root_names.append(root_name)

    # get any text files if any

    for one in txt_files:
        shutil.copy(one,out_dir)
    			

    print('The commands that will be executed will be:')
    for one in commands:
        print(one)

    # Switch to the work directory
    cwd=os.getcwd()
    os.chdir(out_dir)

    proc=subprocess.Popen('Setup_Py_Dir',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    f=open(outputfile,'w')
    f.close()


    run_cmds(commands,root_names,outputfile)

    # Special tests (those which cannot be run as above should be
    # placed here.  The special tests should generally be self
    # contained, that is to say one_liners here.
	
    #  Run a special test on hydro

    # htest=0
    # if os.path.exists(pf_dir+'/hydro'):
    #     print("Hydro directory exists")
    #     hydro_dir=pf_dir+'/hydro'
    #     hydro_files=glob(hydro_dir+'/*')
		
    #     htest=1
		
    # if htest:
    #     for one in hydro_files:
    #         shutil.copy(one,'../'+out_dir)
    #     py_hydro(version,outputfile)

    py_hydro(version, pf_dir, outputfile)

    # run a special macro-atom test -- XXX needs fixing 
    # balmer_decrement.BalmerTest(glob(pf_dir+'/matom_balmer/balmer_test.pf', plotit=True)


    # Return to the place where the code was made from
    os.chdir(cwd)
    return out_dir
	
def run_cmds(commands,root_names,outputfile):
    f=open(outputfile,'a')
    
    i=0
    while i<len(commands):
        one=commands[i]
        string='\n\nRunning %s' % one
        print(string)
        f.write('%s\n'% string)

        proc=subprocess.Popen(one,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr=proc.communicate()
        if len(stderr):
            string='Stderrs were reported!!'
            f.write('%s\n'% string)
            f.write('%s\n' % stderr.decode())
            print(string)
            print(stderr.decode())
            g=open(root_names[i]+'.stderr.txt','w')
            g.write('%s\n' % stderr.decode())
            g.close()
        else:
            string='No stderrs were reported'
            print(string)
            f.write('%s\n' % string)
            # time.sleep(5.)   #  Pause to make sure all of the output files have been written
            check_one(f,root_names[i])

        string='Finished %s' % one
        print(string)
        f.write('%s\n\n'% string)
        f.flush()
        i+=1
    f.close()
		
    return
	
	
	
def py_hydro(version,pf_dir,outputfile):
	
    '''
    This is a self contained script to perform specialised tests of
    the py_hydro models. We need to run on one processor, otherwise
    the test of the heatcool wont work - but its quick!

    Notes:

    This routine is run from within the current working directory,
    which is the directory created earlier be doit


    History

    1801  nsh
    1802  ksl Modified this routine so that it is more standalone
              than previously.  Special regression tests need to
              be isolated from the internal logic of doit so they
              can be added/removed easily.
    2908 kdl  Eliminated the checks for a difference in heating and
              and cooling rates.  If this is important, a check of
              these should be done in regression_checks, and the 
              check should be made between this run and a previous
              run, not one against a file made in the distant past.
    '''

    out_dir=os.getcwd()

    htest=0
    hydro_dir=pf_dir+'/hydro'
    if os.path.exists(hydro_dir):
        print("\nHydro directory exists")
        hydro_files=glob(hydro_dir+'/*')
        htest=1
    else:
        print('The Hydro directory %s doues not appear to exist' % hydro_dir)
        return
		
    if htest:
        for one in hydro_files:
            shutil.copy(one,out_dir)

    root_name=['py_hydro']
    hydro_command=['%s %s  >%s.stdout.txt' % (version,root_name[0]+'.pf',root_name[0])]
    run_cmds(hydro_command,root_name,outputfile)
	
    cmd='cp py_hydro.wind_save py_hydro_restart.wind_save'
    subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    cmd='cp py_hydro.spec_save py_hydro_restart.spec_save'
    subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	
#OLD    cmd='diff py_heatcool.dat model_heatcool.dat'
#OLD    proc=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#OLD    stdout,stderr=proc.communicate()
	
#OLD    if len(stdout):
#OLD        string="py_heatcool.dat has changed from model - important to investigate"
#OLD    else:
#OLD        string="py_heatcool.dat is unchanged - test passed"
#OLD    print(string)
#OLD    f1=open('Summary.txt','a')
#OLD    f1.write('%s\n'% string)
#OLD    f1.close()
	
		
    root_name=['py_hydro_restart']
    # hydro_command=['%s -z -r %s  >%s.stdout.txt' % (version,root_name[0]+'.pf',root_name[0])]
    hydro_command=['%s -r %s  >%s.stdout.txt' % (version,root_name[0]+'.pf',root_name[0])]
    run_cmds(hydro_command,root_name,outputfile)
	
#OLD    cmd='diff py_heatcool.dat model_restart_heatcool.dat'
#OLD    proc=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#OLD    stdout,stderr=proc.communicate()
	
#OLD    if len(stdout):
#OLD        string="restart py_heatcool.dat has changed from model - important to investigate"
#OLD    else:
#OLD        string="restart py_heatcool.dat is unchanged - test passed"
#OLD    print(string)
#OLD    f1=open('Summary.txt','a')
#OLD    f1.write('%s\n'% string)
#OLD    f1.close()


    return



def steer(argv):
    '''
    This is just a steering routine so that switches can be processed 
    from the command line
    '''
    pf_dir=''
    out_dir=''
    np=3

    i=1
    words=[]
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-np':
            i=i+1
            np=int(argv[i])
        elif argv[i]=='-pf_dir':
            i=i+1
            pf_dir=(argv[i])
        elif argv[i]=='-out_dir':
            i=i+1
            out_dir=(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknown switch ---  %s' % argv[i])
            return
        else:
            words.append(argv[i])
        i+=1

    if(len(words)==0):
        print('Error: Consumed of command line without a python executable')
        return

    for one in words:
        q=doit(version=one,pf_dir=pf_dir,out_dir=out_dir,np=np,outputfile='Summary.txt')

    # Now run regression checks between this run and the last time the routine was run

    regression_check.doit(q)







# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print(__doc__)
