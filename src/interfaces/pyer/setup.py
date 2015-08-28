#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension

setup(name="PackageName",
	version='0.1',
	description='calculates the likelihood of relatedness data',
	author='Matthew Ackerman',
	author_email='matthew.s.ackerman@gmail.com',
	url='https://',
	ext_modules=[
		Extension("rml", ["rml/rml.cpp"],
		libraries=["boost_python"],
		extra_compile_args=['-std=c++11', '-fopenmp'],
		extra_link_args=['-lgomp'])
	])
