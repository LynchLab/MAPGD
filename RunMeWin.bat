@echo off
set /p infile=please enter the name of the mpileup file to be analyzed:

echo %infile%

.\bin\mapgd.exe proview -c 5 %infile% > infile.txt
.\bin\mapgd.exe ep -i infile.txt

@echo off
set /p 1=please any key to exit:
