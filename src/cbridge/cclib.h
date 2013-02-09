// Header for C functions to convert python to c

#include <Python.h>

PyObject* getModuleByName(const char* name);
PyObject* getParserModule();
PyObject* getMethodModule();
PyObject* ccopen(PyObject* parseModule, const char* filename);

//eof
