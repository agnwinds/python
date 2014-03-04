#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Check the list of files created for the disk grid to see that they all
exist and make a list of those that are extra.


Command line usage (if any):

	usage: check.py filename

Description:  

Primary routines:

Notes:
									   
History:

140221 ksl Coding begun

'''

import sys
import glob

def read_file(filename,char=''):
	'''
	Read a file and split it into words, eliminating comments
	
	char is an optional parameter used as the delimiter for
	splitting lines into words.  Otherwise white space is
	assumed.

	History:
	
	110729	ksl	Added optional delimiters
	'''

	try:
		f=open(filename,'r')
		xlines=f.readlines()
		f.close()
	except IOError :
		print "The file %s does not exist" % filename
		return []   
	
	lines=[]
	
	i=0
	while i<len(xlines):
		z=xlines[i].strip()
		if char=='':
			z=z.split()
		else:
			z=z.split(char)
		if len(z)>0:
			if z[0][0]!='#':
				lines=lines+[z]
		i=i+1
	return lines


def doit(filename='python.ls'):
	'''
	Do something useful
	'''
	lines=read_file(filename)
	files=glob.glob('*.11')

	print len(lines),len(files)
	python_list=[]
	for line in lines:
		word=line[0].split('/')
		python_list.append(word[1])

	ngood=nbad=0

	print files[0],python_list[0]

	for one in python_list:
		ok='no'
		for  file in files:
			# print 'x',one,file
			if one==file:
			# if one.count(file)>0:
				ok='yes'
				print one,file
				break
		if ok=='yes':
			ngood=ngood+1
		else:
			nbad=nbad+1
			print one
	print 'ngood,nbad',ngood,nbad

	g=open('extra.txt','w')
	extra=[]
	for file in files:
		ok='no'
		for one in python_list:
			if one.count(file)>0:
				ok='yes'
				break
		if ok=='no':
			extra.append(file)
			g.write('rm %s\n' % file)
	g.close()
	


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		# doit(int(sys.argv[1]))
		doit(sys.argv[1])
	else:
		print 'usage: check.py filename'
