#include "CommandLine.hpp"

#include <cstdlib>
#include <cstring>

int main(int argc, char** argv)
{
	atexit([](){std::cout << "Press Enter to exit!"; std::cin.get(); });
	CommandLine cmdl;
	cmdl.run();
	return 0;
}
