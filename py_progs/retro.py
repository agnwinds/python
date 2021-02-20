#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

The functions here are intended to ease running a collection 
of commits on the same parameter file so that one can locate
where a change in Python occured. 



Command line usage (if any):

    None at present. These are intend to be run as part 
    as a Jupyter script or from a python window

Description:  

Primary routines:
    log2table  - create an ascii table that summarizes the log and contains
        entries for compiling and running python many times
    compile_many - create at lot of Python executables in a local directory
    run_many  - runs a set of python executables
    plot_many  - plot consecutive sets of the outputs 


Notes:

    To use this, create a directory where one wants to run

    Create a text file by running 

    git log > commit.txt (or whatever

    Edit this file so it includes only the commit range in which you are interested.

    Run the function log2table to create a file which contains a summary of the logfile.
    This file, hereafter the masterfile, has additional columnns needed for 
    the remainder of the effort

    Run the function compile_many to compile every nth commit.  The various versions of 
    python will be called py_x001, etc, where 001 corresponds to the nubmer of the commit
    going from most recent to oldest.

    Put the .pf file in your local directory.  The runs are going to occur in a subdirectories
    name r_py_x001 etc.  

    If the run requires special files, such as seds that are not part of the normal python 
    distribution, one needs to adjust the parameter files to work in the subdirectories.

    Run the function run_many to run all of the model that you have compiled.  

    Run the function plot_many to examine what the various spectra look like

                                       
History:

210220 ksl Coding begun

'''

import sys
import os
import subprocess
import shutil


from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt



from scipy.signal import boxcar
from scipy.signal import convolve


# Define the directory from which the
# program is called.
home=os.getcwd()


def read_table(filename='foo.txt',format=''):
    '''
    Read a file using astropy.io.ascii and
    return this

    Description:

    Notes:

    This version handles changing string
    lengths better than a simply reading
    the table with ascii.read

    History:


    '''
    try:
        if format=='':
            data=ascii.read(filename)
        else:
            data=ascii.read(filename,format=format)
        for col in data.itercols():
            if col.dtype.kind in 'SU':
                data.replace_column(col.name,col.astype('object'))
    except IOError:
        print ('Error: file %s does not appear to exist' % filename)
        return []

    return data


def log2table(logfile='commits.txt',masterfile='retro_master.txt'):
    '''
    Create a ascii table from a file that contains 
    a normal git log (or the portion of a log that 
    has been edited to be in the time period of interest

    This is normall produced by 

        git log > commits.txt

    One should ensure that the edited file contains complete
    records.
    '''
    x=open(logfile)
    line=x.readlines()
    commit=[]
    author=[]
    date=[]
    info=[]
    xinfo=''
    i=0
    while i < len(line):
        q=line[i].strip()
        word=q.split()
        if len(q)==0:
            pass
        elif word[0]=='commit':
            if len(xinfo):
                if len(xinfo)>100:
                    info.append(xinfo[0:100])
                else:
                    info.append(xinfo)
            xinfo=''
            commit.append(word[1])
            xnext='author'
        elif word[0]=='Author:' and xnext=='author':
            author.append(word[1])
            xnext='date'
        elif word[0]=='Date:' and xnext=='date':
            xdate='%s %s %s %s' %(word[2],word[3],word[4],word[5])
            date.append(xdate)
            xnext='info'
        elif xnext=='info':
            xinfo='%s %s' % (xinfo,q)
        i+=1
    if len(xinfo)>100:
        info.append(xinfo[0:100])
    else:
        info.append(xinfo)
    
    print(len(commit),len(author),len(date),len(info))
    
    xnum=np.linspace(1,len(commit),len(commit),dtype=int)
    
    result=Table([xnum,commit,author,date,info], names=['Number','Commit','Authur','Date','Info'])

    result['Compiled']='No'
    result['Changed']='Unknown'
    result['Run']='No'
    
    result.write(masterfile,format='ascii.fixed_width_two_line',overwrite=True)
        


def compile_one(commit='a77ec180c017244cb56f41b50178daf81541748a',number=35,print_output=True):
    '''
    Compile a single commit of Python and move the executable to a local directory.
    The name of the executable is determined  by the number
    '''
    cwd=os.getcwd()
    os.chdir('/Users/long/Python/source/')
    result=subprocess.run(['git', 'checkout', commit],capture_output=True,text=True)
    if print_output:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)
    result=subprocess.run(['make', 'clean'],capture_output=True,text=True)
    if print_output:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)
    version='VERSION=_x%03d' % number
    
    result=subprocess.run(['make', version, 'python'],capture_output=True,text=True)
    if print_output:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)
    

    if os.path.isfile(cwd+'/py_x%03d' % number):
        os.remove(cwd+'/py_x%03d' % number)
    
    success=True
    try:
        shutil.move('../bin/py_x%03d' % number,cwd)
    except:
        print('Could not find ../bin/py_x%03d for commit %s' %(number,commit))
        success=False

    # return the repository to a known state
    result=subprocess.run(['git', 'checkout', 'dev'],capture_output=True,text=True)

    os.chdir(cwd)
    return success



def check_inputs(command='py',parameter_file='cv.pf',print_output=True):
    '''
    Check that a parameter file has the appropriate parameters to run with 
    a version of python
    '''
    result=subprocess.run([command, '-i', parameter_file],capture_output=True,text=True,timeout=5)
    if print_output:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)
    if result.stdout.count('Exiting'):
        print('Parameter File problem with %s %s' % (command,parameter_file))
        return False
    return True



def check_completion(parameter_file='cv.pf'):
    '''
    Check the .sig file to see if a run has completed
    '''
    root=parameter_file.replace('.pf','')
    xsig=root+'.sig'
    f=open(xsig)
    lines=f.readlines()
    print(lines[-1])
    if lines[-1].count('COMPLETE'):
        return True
    return('False')

def run_one(command='py',parameter_file='cv.pf',np=8, print_output=False,timeout=3000):
    '''
    Run a single model, with a timeout

    Here the timeout does not use the timeout for Python itself but is on the
    python subprocess call.
    '''
    
    if check_inputs(command,parameter_file,print_output)==False:
        return False
    
    try:
        result=subprocess.run(['Setup_Py_Dir'],capture_output=True,text=True)
        result=subprocess.run(['mpirun', '-np',str(np),command, parameter_file],capture_output=True,text=True,timeout=timeout)
    except subprocess.TimeoutExpired:
        return False
    if print_output:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)
    if result.stdout.count('Exiting'):
        print('Parameter File problem with %s %s' % (command,parameter_file))
        return False
    return check_completion(parameter_file)


def compile_many(masterfile='retro_master.txt',increment=10):
    '''
    This routine compiles a number of versions of python and checks
    to see if the code has changed between the time it has changed.
    It adds a flag about this to the table

    '''
    xsum=read_table(masterfile)
    i=0
    old=-1
    while i<len(xsum):
        print(xsum['Commit'][i],xsum['Number'][i])
        success=compile_one(xsum['Commit'][i],xsum['Number'][i],print_output=False)
        if success:
            xsum['Compiled'][i]='Yes'
            new=i
            if old== -1:
                xsum['Changed'][i]='Yes'
            else:
                if os.path.getsize('py_x%03d' % xsum['Number'][old])==os.path.getsize('py_x%03d' % xsum['Number'][i]):
                    xsum['Changed'][i]='No'
                    print('No Change for py_x%03d' % (xsum['Number'][i]) )
                else:
                    xsum['Changed'][i]='Yes'
            old=i
        else:
            xsum['Compiled'][i]='Failed'
            xsum['Changed'][i]='Failed'
        xsum.write(masterfile,format='ascii.fixed_width_two_line',overwrite=True)
        i+=increment
    return 


def run_many(master='retro_master.txt',parameter_file='cv.pf',np=8,rerun=False):
    '''
    Run many versions of python on a single parameter file.
    
    The routines are run in a subdirectory that is defined by the
    version number, so the parameter files must be able to run from there.
    In particular if any special files are called then the path to the
    files must be defined correctly, and any set of model files must
    be defined so that they will be found from the subdirectory as well.
    '''
    start_dir=os.getcwd()
    x=read_table(master)
    # Check that the file contains the column Ran, and if not add it
    check=False
    for one in x.colnames:
        if one=='Ran':
            check=True
    if rerun or check==False:
        x['Ran']='           No'
    
    for one in x:
        if one['Changed']=='Yes' and one['Ran']!='Yes':
            xcommand='py_x%03d' % one['Number']
            if os.path.isdir('r_'+xcommand)==False:
                os.mkdir('r_'+xcommand)
            shutil.copy(parameter_file,'r_'+xcommand)
            os.chdir('r_'+xcommand)
            success=run_one(command='../'+xcommand,parameter_file=parameter_file,np=np, print_output=False,timeout=8000)
            if success==True:
                one['Ran']='Yes'
            else:
                one['Ran']='Failed'
            os.chdir(start_dir)
            x.write(master,format='ascii.fixed_width_two_line',overwrite=True)
    return 



def xsmooth(flux,smooth=21):
    '''
    boxcar smooth the flux

    This is just a little utility to 
    allow me  to boxcar smoothe thins
    without having to call the appropirate
    routines from scipy
    '''
    if (smooth)>1:
        q=convolve(flux,boxcar(smooth)/float(smooth),mode='same')
        return(q)
    else:
        return(flux)



def plot_two(dir1='r_py_x001',dir2='r_py_x026',file='agn.spec',colname='A01P0.50',smooth=1,xlim=[0,0],ylim=[0,0]):
    '''
    Create plots that compare two Python runs
    '''
    
    try:
        one=ascii.read('%s/%s' %(dir1,file))
    except:
        print('Could not read %s/%s' %(dir1,file))
        return 1
    try:
        two=ascii.read('%s/%s' %(dir2,file))
    except:
        print('Could not read %s/%s' %(dir2,file))   
        return 2
    
    got_col=False
    for col in one.colnames:
        if colname==col:
            got_col=True
    if got_col==False:
        print('File does not have col %s' % colname)
        print('The possible column names are :\n',one.colnames)
        return 1

    root1=dir1.replace('r_py_','')
    root2=dir2.replace('r_py_','')

    if smooth>1:
        y1=xsmooth(one[colname],smooth)
        y2=xsmooth(two[colname],smooth)
    else:
        y1=one[colname]
        y2=two[colname]


    
    plt.figure(1)
    plt.clf()
    plt.plot(one['Lambda'],y1,label=root1)
    plt.plot(two['Lambda'],y2,label=root2)
    if xlim[0]!=xlim[1]:
        plt.xlim(xlim[0],xlim[1])
    if ylim[0]!=ylim[1]:
        plt.ylim(ylim[0],ylim[1])
    plt.legend()
    
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.tight_layout()
    plt.savefig('%s_%s.png' % (root1,root2))
    return 0






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print ('usage: xsmooth filename')

def plot_many(masterfile='retro_master.txt',file='agn.spec',colname='A01P0.50',smooth=11,xlim=[0,0],ylim=[0,0]):
    '''
    Plot spectra from a set of models all of which should h
    '''
    x=ascii.read(masterfile)
    xx=x[x['Ran']=='Yes']
    print(len(xx))
    i=1
    while i<len(xx):
        dir1='r_py_x%03d' % (xx['Number'][i-1])
        dir2='r_py_x%03d' % (xx['Number'][i])
        plot_two(dir1,dir2,file,colname,smooth,xlim,ylim)
        i+=1
    return
        

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print ('usage: retro.py filename')
