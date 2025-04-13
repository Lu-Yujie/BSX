#ifndef SUBGRAPHMATCHING_MATCHINGCOMMAND_H
#define SUBGRAPHMATCHING_MATCHINGCOMMAND_H

#include "utility/commandparser.h"
#include <map>
#include <iostream>
enum OptionKeyword {
    QueryGraphFile,             // -q, The query graph file path, compulsive parameter
    DataGraphFile,              // -d, The data graph file path, compulsive parameter
    Filter,                     // -filter, The strategy of filtering
    Order,                      // -order, The strategy of ordering
    Engine,                     // -engine, The computation engine
    MaxOutputEmbeddingNum,      // -num, The maximum output embedding num
    TimeLimit,                  // -time_limit, The time limit for executing a query in milliseconds
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

    std::string getFilterType() {
        return options_value[OptionKeyword::Filter] == "" ? "CFL" : options_value[OptionKeyword::Filter];
    }

    std::string getOrderType() {
        return options_value[OptionKeyword::Order] == "" ? "GQL" : options_value[OptionKeyword::Order];
    }

    std::string getEngineType() {
        return options_value[OptionKeyword::Engine] == "" ? "LFTJ" : options_value[OptionKeyword::Engine];
    }

    std::string getMaximumEmbeddingNum() {
        return options_value[OptionKeyword::MaxOutputEmbeddingNum] == "" ? "MAX" : options_value[OptionKeyword::MaxOutputEmbeddingNum];
    }

    std::string getTimeLimit() {
        return options_value[OptionKeyword::TimeLimit] == "" ? "1000" : options_value[OptionKeyword::TimeLimit];
    }

};


#endif //SUBGRAPHMATCHING_MATCHINGCOMMAND_H
