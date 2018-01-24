#!/usr/bin/env python 


'''
                    Space Telescope Science Institute

Synopsis:  

Replace current informational headers with one that can
be parsed with doxygen


Command line usage (if any):

    usage: doxygen.py whatever.c

Description:  

    This is an attempt to write a parser which will read one of the
    .c files used in Python and locate the headers that are supposed
    to be at the top of every function call.  Parsing these headers
    and also using information from cproto the routine attemps to 
    create a partially populated doxygen block for each routine.

    The routine writes out a new file, new_whatever.c containging
    doxygen blocks just in fromt of every function call.

    It is up to the user to copy new_whatever.c back to whatever.c once
    he/she has thoroughly looked at the file.

Primary routines:

    doit

Notes:

    The routine has to be run in a directory that contains the header 
    files, e.g atomic.h, so that cproto will perform properly.

    The routine will not work with 'non-standard' old headers, so
    one should carefully look at all of the created blocks to see that
    they make sense.

    The old headers are commented out with //OLD and at once we are fairly
    sure they are no longer needed.

    If there was not old_header, there will still be a partially filled out
    doxygen block in front of each routine in the file

    Right now there is a lot of repetive code for populating the various 
    blocks.  ksl suspects this should be systemetized


                                       
History:

180123 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy
import subprocess


def read_file(filename,char=''):
    '''
    Read a file
    
    '''

    try:
        f=open(filename,'r')
        xlines=f.readlines()
        f.close()
    except IOError :
        print ("The file %s does not exist" % filename)
        return []   

    return xlines


def get_modules(filename='emission.c'):
    '''
    use cproto to capture the functions etc that 
    are contained in the routine.  Split what is
    returned into something which can be incorporated
    into a search

    '''
    proc=subprocess.Popen('cproto %s ' % filename,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()

    stdout=stdout.decode()
    stdout=stdout.split('\n')

    print(stdout)

    records=[]
    i=1
    while i<len(stdout)-1:
        z=stdout[i].replace('(',' ')
        z=z.replace(')',' ')
        z=z.replace(',',' ')
        z=z.replace(';',' ')
        words=z.split()
        print(words)
        one_record=[stdout[i]]+words
        records.append(one_record)
        # records.append([stdout[i],words[0],words[1]])

        i+=1

    return(records)

file_string='''

/***********************************************************/
/** @file   reverb.c
* @Author SWM
* @date   July, 2015
* @brief  Reverberation mapping functions.
*
* File containing reverberation mapping functions.
***********************************************************/

'''

module_string_start='''

/**********************************************************/
/** @name       %s
* @brief       %s
*
'''
module_string_param='''* @param [in]  %s  %s   ???
'''
module_string_end='''
* @return       %s
*
* %s
*
* @notes
* 
* %s
*
***********************************************************/

'''

def doit(filename='emission.c',outputfile=''):
    '''
    Do something magnificent

    Description:

    Notes:

    For cproto to work properly, files like atomic.h and python.h
    must be in the same directory

    History:


    '''

    if outputfile=='':
        outputfile='new_'+filename

    lines=read_file(filename)

    if len(lines)==0:
        return

    modules=get_modules(filename)

    print(modules)



    xlines=0
    i=0
    j=0
    header_start=[]
    header_end=[]
    module_start=[]
    while i<len(lines)-1:
        line=lines[i]
        # print(line)
        header=False

        if line[0:24]=='/***********************':
            print('This is the beginning of a header')
            header_start.append(i)
            header=True
        if line[0:20]=='********************':
            print('Found the end of a header')
            header_end.append(i)
            header=False
        if header==False and j<len(modules):
            # Look for the begining of a module
            x=line.strip()
            y=lines[i+1].strip()
            print('test',x,y)
            x=x.split()
            y=y.split()
            print('test2',x,'y',y)
            if len(x)>0 and len(y)>0 and x[0]==modules[j][1] and y[0]==modules[j][2]:
                module_start.append(i)
                print('Found module start')
                j+=1

        i+=1

    print(header_start)
    print(header_end)
    print(module_start)


    # Now we need to try to match the current headers with the modules because
    # some may be missing

    xmatch=[]
    i=0
    while i<len(modules):
        xmatch.append(-1)
        i+=1

    j=0
    isep=20
    while j<len(header_end):
        k=0
        while k<len(module_start):
            if header_end[j]<module_start[k] and module_start[k]-header_end[j]<isep:
                xmatch[k]=j
                break
            k+=1
        j+=1

    print('xmatch',xmatch)


    x=open(outputfile,'w')

    i=0
    kk=0
    while i<len(lines):
        line=lines[i]
        k=0
        while k<len(header_end):
            if header_start[k]<=i and i<= header_end[k]:
                line='//OLD '+line
                kcurrent=k
                break
            k+=1

        if kk<len(module_start) and i==module_start[kk]:
            proto_string=modules[kk]
            if xmatch[kk]!=-1:
                istart=header_start[kcurrent]
                istop=header_end[kcurrent]
                # Try to get the synopusis
                synopsis_string=''
                synopsis=False
                ii=istart
                while ii<istop:
                    if lines[ii].count('Synopsis:'):
                        print('gotcha')
                        synopsis_string=lines[ii].replace('Synopsis:','')
                        synopsis=True
                    elif lines[ii].count('Arguments:') or lines[ii].count('Returns:'):
                        synopsis=False
                    elif synopsis==True:
                        synopsis_string=synopsis_string+lines[ii]
                    ii+=1
            print('test',synopsis_string)

            if synopsis_string=='':
                synopsis_string='???'
            else:
                synopsis_string=synopsis_string.replace('\n','\n*      ')
            x.write(module_string_start % (proto_string[2],synopsis_string))


        # Now try to get the return string
        if kk<len(module_start) and i==module_start[kk]:
            proto_string=modules[kk]
            if xmatch[kk]!=-1:
                istart=header_start[kcurrent]
                istop=header_end[kcurrent]
                # Try to get the synopusis
                return_string=''
                xreturn=False
                ii=istart
                while ii<istop:
                    if lines[ii].count('Returns:'):
                        print('gotcha')
                        return_string=lines[ii].replace('Returns:','')
                        xreturn=True
                    elif lines[ii].count('Description:'):
                        xreturn=False
                    elif xreturn==True:
                        return_string=return_string+lines[ii]
                    ii+=1
            print('test',return_string)

            return_string=return_string.strip()
            if return_string=='':
                return_string='???'
            else:
                return_string=return_string.replace('\n','\n*      ')


        # Now try to get the description   
        if kk<len(module_start) and i==module_start[kk]:
            proto_string=modules[kk]
            if xmatch[kk]!=-1:
                istart=header_start[kcurrent]
                istop=header_end[kcurrent]
                # Try to get the synopusis
                description_string=''
                description=False
                ii=istart
                while ii<istop:
                    if lines[ii].count('Description:'):
                        print('gotcha')
                        description_string=lines[ii].replace('Description:','')
                        description=True
                    elif lines[ii].count('Notes:'):
                        description=False
                    elif description==True:
                        description_string=description_string+lines[ii]
                    ii+=1
            print('test',description_string)

            description_string=description_string.strip()
            if description_string=='':
                description_string='???'
            else:
                description_string=description_string.replace('\n','\n*      ')



        # Now try to get the notes   
        if kk<len(module_start) and i==module_start[kk]:
            proto_string=modules[kk]
            if xmatch[kk]!=-1:
                istart=header_start[kcurrent]
                istop=header_end[kcurrent]
                # Try to get the synopusis
                notes_string=''
                notes=False
                ii=istart
                while ii<istop:
                    if lines[ii].count('Notes:'):
                        print('gotcha')
                        notes_string=lines[ii].replace('Notes:','')
                        notes=True
                    elif lines[ii].count('History:'):
                        notes=False
                    elif notes==True:
                        notes_string=notes_string+lines[ii]
                    ii+=1
            print('Notes_test',notes_string)

            notes_string=notes_string.strip()
            if notes_string=='':
                notes_string='None'
            else:
                notes_string=notes_string.replace('\n','\n*      ')

                

            kkk=2
            print(proto_string)
            while kkk<len(proto_string):
                x.write(module_string_param % (proto_string[kkk-1],proto_string[kkk]))
                kkk+=2
            x.write(module_string_end % (return_string,description_string,notes_string))
            kk+=1

        
        x.write(line)
        i+=1






    

    return


    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
     
    if len(sys.argv)==1 or sys.argv[1]=='-h':
        print(__doc__)
    elif len(sys.argv)>1:
        doit(sys.argv[1])
    else:
        print ('usage: doxygen.py [-h] filename')
