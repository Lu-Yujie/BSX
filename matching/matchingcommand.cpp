//
// Created by Shixuan Sun on 2018/6/29.
//

#include "matchingcommand.h"

MatchingCommand::MatchingCommand(const int argc, char **argv) : CommandParser(argc, argv) {
    // Initialize options value
    options_key[OptionKeyword::Algorithm] = "-a";
    options_key[OptionKeyword::IndexType] = "-i";
    options_key[OptionKeyword::QueryGraphFile] = "-q";
    options_key[OptionKeyword::DataGraphFile] = "-d";
    options_key[OptionKeyword::ThreadCount] = "-n";
    options_key[OptionKeyword::DepthThreshold] = "-d0";
    options_key[OptionKeyword::WidthThreshold] = "-w0";
    options_key[OptionKeyword::Filter] = "-filter";
    options_key[OptionKeyword::Order] = "-order";
    options_key[OptionKeyword::Engine] = "-engine";
    options_key[OptionKeyword::MaxOutputEmbeddingNum] = "-num";
    options_key[OptionKeyword::ExeTimeLimit] = "-time_limit";
    options_key[OptionKeyword::CSRFilePath] = "-csr";
    options_key[OptionKeyword::OutputFile] = "-o";
    options_key[OptionKeyword::ConfPath] = "-conf";
    processOptions();
};

void MatchingCommand::processOptions() {
    // Query graph file path
    options_value[OptionKeyword::QueryGraphFile] = getCommandOption(options_key[OptionKeyword::QueryGraphFile]);;

    // Data graph file path
    options_value[OptionKeyword::DataGraphFile] = getCommandOption(options_key[OptionKeyword::DataGraphFile]);

    // Algorithm
    options_value[OptionKeyword::Algorithm] = getCommandOption(options_key[OptionKeyword::Algorithm]);

    // Thread count
    options_value[OptionKeyword::ThreadCount] = getCommandOption(options_key[OptionKeyword::ThreadCount]);

    // Depth threshold
    options_value[OptionKeyword::DepthThreshold] = getCommandOption(options_key[OptionKeyword::DepthThreshold]);

    // Width threshold
    options_value[OptionKeyword::WidthThreshold] = getCommandOption(options_key[OptionKeyword::WidthThreshold]);

    // Index Type
    options_value[OptionKeyword::IndexType] = getCommandOption(options_key[OptionKeyword::IndexType]);

    // Filter Type
    options_value[OptionKeyword::Filter] = getCommandOption(options_key[OptionKeyword::Filter]);

    // Order Type
    options_value[OptionKeyword::Order] = getCommandOption(options_key[OptionKeyword::Order]);

    // Engine Type
    options_value[OptionKeyword::Engine] = getCommandOption(options_key[OptionKeyword::Engine]);

    // Maximum output embedding num.
    options_value[OptionKeyword::MaxOutputEmbeddingNum] = getCommandOption(options_key[OptionKeyword::MaxOutputEmbeddingNum]);

    // Time Limit
    options_value[OptionKeyword::ExeTimeLimit] = getCommandOption(options_key[OptionKeyword::ExeTimeLimit]);

    // CSR file path
    options_value[OptionKeyword::CSRFilePath] = getCommandOption(options_key[OptionKeyword::CSRFilePath]);

    // output file path
    options_value[OptionKeyword::OutputFile] = getCommandOption(options_key[OptionKeyword::OutputFile]);

    // configuration files path
    options_value[OptionKeyword::ConfPath] = getCommandOption(options_key[OptionKeyword::ConfPath]);
}