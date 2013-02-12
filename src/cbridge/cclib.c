// Implementations for cclib bridge

#include <Python.h>
#include "cclib.h"

PyObject* ccopen(const char* filename) {

    PyGILState_STATE gstate;
    gstate = PyGILState_Ensure();

    PyObject *parserModule, *ccopenFunction, *filenameObject, *parserObject = NULL;
    PyObject *args;

    if (!filename) {
        printf("Invalid filename.\n");
        return NULL;
    }

    parserModule = PyImport_ImportModule("cclib.parser");
    if ( !parserModule ) {
        printf("cclib.parser not imported.\n");
        PyGILState_Release(gstate);
        return NULL;

    }

    ccopenFunction = PyObject_GetAttrString(parserModule, "ccopen");
    if ( ccopenFunction && PyCallable_Check(ccopenFunction)) {
    	filenameObject = PyString_FromString(filename);
        if (!filenameObject) {
            printf("There's a problem with %s.\n", filename);
            Py_DECREF(ccopenFunction);
            Py_DECREF(parserModule);
            PyGILState_Release(gstate);
            return NULL;
        }

        args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, filenameObject);
        parserObject = PyObject_CallObject(ccopenFunction, args);
        Py_DECREF(args);

    }

    else {
        printf("There is a problem with ccopen.\n");

    }

    PyGILState_Release(gstate);
    printf("Returning parser object.\n");
    return parserObject;
}

