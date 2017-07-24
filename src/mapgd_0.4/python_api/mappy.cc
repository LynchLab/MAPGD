#ifndef MAPGDAPI
#define MAPGDAPI
#include <Python.h>
#include <iostream>

//#include "pymap.h"

#include "map-file.h"
#include "data.h"

struct file_s{
	Base_file *file;
	Data *data;
	std::fstream *fs;
	void kill_me (void *)
	{
		delete file;
		delete data;
		delete fs;
	};
};

static PyObject * 
open(PyObject *self, PyObject *args, PyObject *kwargs)
{
	const char *filename;
	int mode = 2;
	static char *keywords[] = {"filename", "mode", NULL};
	PyArg_ParseTuple(args, "s", &filename);

	file_s *f=new file_s;
	f->fs=new std::fstream;
	f->fs->open(filename);
	f->file=new Base_file;
	f->file->open((std::istream *)(f->fs), READ);
        f->data=f->file->read_header();
	
	return PyCObject_FromVoidPtr( (void*) f, (void (*)(void*))(&file_s::kill_me)  );
}

static PyObject *
read_line(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyObject *obj;
	PyArg_ParseTuple(args, "O", &obj);
	file_s *f=(file_s *) PyCObject_AsVoidPtr( obj );
	std::stringstream out;
	if ( f->file->is_open() )
	{
		while( f->file->read(f->data).table_is_open() ) out << *f->data << std::endl;
	}
	return Py_BuildValue("s", out.str().c_str() );
	std::cerr << "File is not open\n";
	return NULL;
}


static PyMethodDef
module_functions[] = {
	{ "open", (PyCFunction)open, METH_VARARGS | METH_KEYWORDS, "open" },
	{ "read", (PyCFunction)read_line, METH_VARARGS | METH_KEYWORDS, "read" },
	{ NULL }
};

/*
PyMODINIT_FUNC
PyInit_mapgd(void)
{
	PyObject *m;
	static void *PySpam_API[PySpam_API_pointers];
	PyObject *c_api_object;

	m = PyModule_Create(&spammodule);
	if (m == NULL)
	return NULL;

	PySpam_API[PySpam_System_NUM] = (void *)PySpam_System;

	c_api_object = PyCapsule_New((void *)PySpam_API, "spam._C_API", NULL);

	if (c_api_object != NULL)
	PyModule_AddObject(m, "_C_API", c_api_object);
	return m;
}*/


PyMODINIT_FUNC
initmappy(void)
{
	Py_InitModule("mappy", module_functions);
}

/*int 
main (int argc, char **argv)
{
	return 0;
};*/
#endif
