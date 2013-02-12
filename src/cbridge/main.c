#include <Python.h>
#include "cclib.h"

int main(int argc, char* argv[]) {

    PyObject *parserObject, *parseFunction, *args, *dataObject;
    
    Py_Initialize();

    parserObject = ccopen(argv[1]);
    if(!parserObject)
        printf("cclib parser object not created.\n");

    parseFunction = PyObject_GetAttrString(parserObject, "parse");
    if( parseFunction && PyCallable_Check(parseFunction)) {
        args = PyTuple_New(0);
        dataObject = PyObject_CallObject(parseFunction, args);
    }

    Py_Finalize();

    printf("Exiting.\n");
    return 0;
}

