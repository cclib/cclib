#include <Python.h>
#include "cclib.h"

int main(int argc, char* argv[]) {

    PyObject *parserModule;
    PyObject *parserObject;
    PyObject *name, *text, *args;

    Py_Initialize();

    parserModule = getParserModule();
    parserObject = ccopen(parserModule, argv[1]);

    name = PyObject_GetAttrString(parserObject, "parse");
    args = PyTuple_New(0);
    text = PyObject_CallObject(name, args);
//    printf( "%s\n", PyString_AsString(text));

    Py_Finalize();
    return 0;
}

