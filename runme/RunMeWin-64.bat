@echo off
set /p infile=please enter the name of the mpileup file to be analyzed:

echo %infile%

.\bin\mapgd-64.exe proview %infile% > infile.txt
.\bin\mapgd-64.exe ep

@echo off
set /p 1=please any key to exit:
