#include "matchingcommand.h"

MatchingCommand::MatchingCommand(const int argc, char **argv) : CommandParser(argc, argv) {
    // Initialize options value
    options_key[OptionKeyword::QueryGraphFile] = "-q";
    options_key[OptionKeyword::DataGraphFile] = "-d";
    options_key[OptionKeyword::Filter] = "-filter";
    options_key[OptionKeyword::Order] = "-order";
    options_key[OptionKeyword::Engine] = "-engine";
    options_key[OptionKeyword::MaxOutputEmbeddingNum] = "-num";
    options_key[OptionKeyword::TimeLimit] = "-time_limit";
    processOptions();
};

void MatchingCommand::processOptions() {
    // Query graph file path
    options_value[OptionKeyword::QueryGraphFile] = getCommandOption(options_key[OptionKeyword::QueryGraphFile]);;

    // Data graph file path
    options_value[OptionKeyword::DataGraphFile] = getCommandOption(options_key[OptionKeyword::DataGraphFile]);

    // Filter Type
    options_value[OptionKeyword::Filter] = getCommandOption(options_key[OptionKeyword::Filter]);

    // Order Type
    options_value[OptionKeyword::Order] = getCommandOption(options_key[OptionKeyword::Order]);

    // Engine Type
    options_value[OptionKeyword::Engine] = getCommandOption(options_key[OptionKeyword::Engine]);

    // Maximum output embedding num.
    options_value[OptionKeyword::MaxOutputEmbeddingNum] = getCommandOption(options_key[OptionKeyword::MaxOutputEmbeddingNum]);

    // Time Limit
    options_value[OptionKeyword::TimeLimit] = getCommandOption(options_key[OptionKeyword::TimeLimit]);
}