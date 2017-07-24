from distutils.core import setup, Extension
import subprocess

batcmd="find ../ -name \"*.o\" -not -name mapgd.o -not -name mappy.o -not -name gzstream.o | tr -s \'\\n\'  \' \' "
result = subprocess.check_output(batcmd, shell=True)
print result
objects=result.split()

minor="a"#open("../VERSION").read()

#larg="-lgsl -lgslcblas -lm -lgzstream -lz -L ./gzstream/"

module1 = Extension('mappy',
		sources = ['mappy.cc'],
		extra_compile_args = ["--std=c++11", "-D VERSION=\"0.4."+minor+"\"", "-fopenmp"],
		extra_link_args=['-lgomp'],
		library_dirs=["../gzstream/", "../../../../htslib"],
		libraries=['gsl', 'gslcblas', "m", "gzstream", "z", "hts", "sqlite3"],
		include_dirs=['../', '../data_types/', '../commands/', '../raw/', '../data_conversion/', '../gzstream/', '../stream_tools/'],
		extra_objects = objects )


setup (name = 'mappy',
	version = "0.4."+minor,
      	author='Matthew Ackerman',
	description = 'This is a mapgd api',
	url='https://www.github.com/LynchLab/mapgd/',
#	packages=['distutils', 'distutils.command'],
	ext_modules = [module1])
