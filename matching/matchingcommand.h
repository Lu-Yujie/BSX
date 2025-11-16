#ifndef SUBGRAPHMATCHING_MATCHINGCOMMAND_H
#define SUBGRAPHMATCHING_MATCHINGCOMMAND_H

#include "commandparser.h"
#include <map>
#include <iostream>
enum OptionKeyword {
    QueryGraphFile,             // -q, The query graph file path, compulsive parameter
    DataGraphFile,              // -d, The data graph file path, compulsive parameter
    MaxOutputEmbeddingNum,      // -num, The maximum output embedding num
    TimeLimit,                  // -time_limit, millisecond
    OutputFile,                 // -o, output file path(absolute)
    ConfPath                    // -conf, configuration file path(absolute) 
    };

class MatchingCommand : public CommandParser{
private:
    std::map<OptionKeyword, std::string> options_key;
    std::map<OptionKeyword, std::string> options_value;

private:
    void processOptions();

public:
    MatchingCommand(int argc, char **argv);

    std::string getDataGraphFilePath() {
        return options_value[OptionKeyword::DataGraphFile];
    }

    std::string getQueryGraphFilePath() {
        return options_value[OptionKeyword::QueryGraphFile];
    }

    std::string getMaximumEmbeddingNum() {
        return options_value[OptionKeyword::MaxOutputEmbeddingNum] == "" ? "MAX" : options_value[OptionKeyword::MaxOutputEmbeddingNum];
    }

    std::string getTimeLimit() {
        return options_value[OptionKeyword::TimeLimit] == "" ? "1" : options_value[OptionKeyword::TimeLimit];
    }

    std::string getOutputFile() {
        return options_value[OptionKeyword::OutputFile] == "" ? "./out_default.csv" : options_value[OptionKeyword::OutputFile];
    }

    std::string getConfPath() {
        return options_value[OptionKeyword::ConfPath] == "" ? "./conf" : options_value[OptionKeyword::ConfPath];
    }
};


#endif //SUBGRAPHMATCHING_MATCHINGCOMMAND_H
