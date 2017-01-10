#include "CommandLine.hpp"

#include "MyApp.h"

#include <cstdlib>
#include <cstring>

int main(int argc, char** argv)
{
	atexit([](){std::cout << "Press Enter to exit!"; std::cin.get(); });
	//if(argc > 1 && strcmp(argv[1], "-gui") == 0)
	{
		// start gui
		CMyApp app;
		app.initGraphics();
		app.Run();
		return 0;
	}
	//else
	{
		CommandLine cmdl;
		cmdl.run();
		return 0;
	}
}
