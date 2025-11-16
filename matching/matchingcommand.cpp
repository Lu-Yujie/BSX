#include "matchingcommand.h"

MatchingCommand::MatchingCommand(const int argc, char **argv) : CommandParser(argc, argv) {
    // Initialize options value
    options_key[OptionKeyword::QueryGraphFile] = "-q";
    options_key[OptionKeyword::DataGraphFile] = "-d";
    options_key[OptionKeyword::MaxOutputEmbeddingNum] = "-num";
    options_key[OptionKeyword::TimeLimit] = "-time_limit";
    options_key[OptionKeyword::OutputFile] = "-o";
    options_key[OptionKeyword::ConfPath] = "-conf";
    processOptions();
};

void MatchingCommand::processOptions() {
    // Query graph file path
    options_value[OptionKeyword::QueryGraphFile] = getCommandOption(options_key[OptionKeyword::QueryGraphFile]);;

    // Data graph file path
    options_value[OptionKeyword::DataGraphFile] = getCommandOption(options_key[OptionKeyword::DataGraphFile]);

    // Maximum output embedding num
    options_value[OptionKeyword::MaxOutputEmbeddingNum] = getCommandOption(options_key[OptionKeyword::MaxOutputEmbeddingNum]);

    // Time Limit
    options_value[OptionKeyword::TimeLimit] = getCommandOption(options_key[OptionKeyword::TimeLimit]);

    // Output file
    options_value[OptionKeyword::OutputFile] = getCommandOption(options_key[OptionKeyword::OutputFile]);

    // Configuartion file path
    options_value[OptionKeyword::ConfPath] = getCommandOption(options_key[OptionKeyword::ConfPath]);
}