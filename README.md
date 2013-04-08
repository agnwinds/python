README
***
=========
Precursor - python_74b5_JM

* This incorporates files which Nick thought were not in python_74b5, but the only one that was really different is the wind2d.c, which is I believe the adiabatic fix

Other small changes are as follows:
* Deleted a number of print statements that Nick had put in for diagnostic purposes.  These are cluttering up the printed outputs which were intended to be brief and make it possible for one to use the printed output to see where one was in the routine.  If needed, one should use the Log or Log_silent so they will go in the diagnostic.  Otherwise ddd is a better approach.  changed others that were Log to Log_silent
* Added error check concerning NDIM_MAX.  The program should exit with an error if this a problem.  It would be possible to define these arrays dynamically, and that is probably what should be done.  In that sense, my fix is a kluge
* Note - This needs regression testing.  Also, it proabably should be chieked to see if the fixes that have been worked out are all incorproated.

