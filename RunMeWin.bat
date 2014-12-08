@echo off
set /p infile=please enter the name of the mpileup file to be analyzed:
./bin/mapgd.exe proview @infile > infile.txt
./bin/mapgd.exe ep -i @infile
