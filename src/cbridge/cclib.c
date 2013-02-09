// Implementations for cclib bridge

#include <Python.h>
#include "cclib.h"

PyObject* getModuleByName(const char* name) {
    PyObject *moduleName, *parserModule;

    moduleName = PyString_FromString(name);
    if( !moduleName )
        return NULL;

    parserModule = PyImport_Import(moduleName);
    Py_DECREF(moduleName);

    return parserModule;
}

PyObject* getParserModule(const char* name) {
    return getModuleByName("cclib.parser");
}

PyObject* getMethodModule(const char* name) {
    return getModuleByName("cclib.method");
}

PyObject* ccopen(PyObject* parserModule, const char* filename) {

    PyObject *ccopenFunction, *filenameObject, *parserObject = NULL;
    PyObject *args;

    if ( !parserModule )
        return NULL;

    ccopenFunction = PyObject_GetAttrString(parserModule, "ccopen");
    if ( ccopenFunction && PyCallable_Check(ccopenFunction)) {
    
        filenameObject = PyString_FromString(filename);
        if (!filenameObject)
            return NULL;

        args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, filenameObject);
        parserObject = PyObject_CallObject(ccopenFunction, args);

    }

    return parserObject;
}

