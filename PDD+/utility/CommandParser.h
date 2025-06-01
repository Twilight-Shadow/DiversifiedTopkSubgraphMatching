#ifndef COMMANDPARSER_H
#define COMMANDPARSER_H

#include <string>
#include <algorithm>
#include <vector>

class CommandParser {
private:
    std::vector<std::string> tokens_;

public:
    CommandParser(const int argc, char **argv);
    const std::string getCommandOption(const std::string &option) const;
    bool commandOptionExists(const std::string &option) const;
};


#endif //COMMANDPARSER_H