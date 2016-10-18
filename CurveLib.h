// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the CURVELIB_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// CURVELIB_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef CURVELIB_EXPORTS
#define CURVELIB_API __declspec(dllexport)
#else
#define CURVELIB_API __declspec(dllimport)
#endif

// This class is exported from the CurveLib.dll
class CURVELIB_API CCurveLib {
public:
	CCurveLib(void);
	// TODO: add your methods here.
};

extern CURVELIB_API int nCurveLib;

CURVELIB_API int fnCurveLib(void);
