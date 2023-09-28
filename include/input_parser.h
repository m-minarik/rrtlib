// https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
// @author iain

#ifndef __INPUT_PARSER__
#define __INPUT_PARSER__

#include <string>
#include <vector>

class InputParser
{
public:
	InputParser(const int &argc, char **argv);
	const std::string &getCmdOption(const std::string &option) const;
	bool cmdOptionExists(const std::string &option) const;

private:
	std::vector<std::string> tokens;
};

#endif // __INPUT_PARSER__
