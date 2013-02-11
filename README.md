README
***
=========
Precursor is python_74b
Changes this version
All sane checks in all files now have an error message associtated with them that reports SUBROUTINE:sane_check <information>
A lot of files have been modified, but those without any errors associated with sane checks were
variable_temperature checks on densities 
thierry.c - few sane checks had errors - but I doubt this file is called very much - seems experimental
vector.c - there was one sane check without an error, for calculating the length of a vector
saha.c - checks on densities - this was the original pattern for the vriable_temperatures code, hence the copied lack of errors!


