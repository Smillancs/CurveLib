// CurveLib.cpp : Defines the exported functions for the DLL application.
//

#include "CurveLib.h"


// This is an example of an exported variable
CURVELIB_API int nCurveLib=0;

// This is an example of an exported function.
CURVELIB_API int fnCurveLib(void)
{
    return 42;
}

// This is the constructor of a class that has been exported.
// see CurveLib.h for the class definition
CCurveLib::CCurveLib()
{
    return;
}
