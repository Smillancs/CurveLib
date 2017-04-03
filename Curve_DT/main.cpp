#include "MyApp.h"

#include <iostream>

int main(int argc, char** argv)
{
	atexit([](){std::cout << "Press Enter to exit!"; std::cin.get(); });
	// start gui
	CMyApp app;
	app.initGraphics();
	app.Run();
	return 0;
}
