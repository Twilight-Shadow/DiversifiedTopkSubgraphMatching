#include "CommandList.h"

MatchingCommand::MatchingCommand(const int argc, char **argv) : CommandParser(argc, argv) {
    // Initialize options value
    options_key[OptionKeyword::QueryGraphFile] = "-q";
    options_key[OptionKeyword::PatitionedDataGraphDirectory] = "-d";
    options_key[OptionKeyword::ThreadCount] = "-n";
    options_key[OptionKeyword::k] = "-k";
    options_key[OptionKeyword::DataGraphFile] = "-g";
    options_key[OptionKeyword::PartitionCount] = "-p";
    options_key[OptionKeyword::DistanceFile] = "-m";
    options_key[OptionKeyword::GroundTruthFile] = "-t";
    options_key[OptionKeyword::ReplicatedNodeDirectory] = "-r";
    processOptions();
};

void MatchingCommand::processOptions() {
    // Query graph file path
    options_value[OptionKeyword::QueryGraphFile] = getCommandOption(options_key[OptionKeyword::QueryGraphFile]);;

    // Partitioned data graph directory path
    options_value[OptionKeyword::PatitionedDataGraphDirectory] = getCommandOption(options_key[OptionKeyword::PatitionedDataGraphDirectory]);

    // Thread count
    options_value[OptionKeyword::ThreadCount] = getCommandOption(options_key[OptionKeyword::ThreadCount]);

    // Top-k value
    options_value[OptionKeyword::k] = getCommandOption(options_key[OptionKeyword::k]);

    // Data graph file path
    options_value[OptionKeyword::DataGraphFile] = getCommandOption(options_key[OptionKeyword::DataGraphFile]);

    // Partition number
    options_value[OptionKeyword::PartitionCount] = getCommandOption(options_key[OptionKeyword::PartitionCount]);

    // Distance file path
    options_value[OptionKeyword::DistanceFile] = getCommandOption(options_key[OptionKeyword::DistanceFile]);

    // Ground truth file path
    options_value[OptionKeyword::GroundTruthFile] = getCommandOption(options_key[OptionKeyword::GroundTruthFile]);

    // Replicated node directory
    options_value[OptionKeyword::ReplicatedNodeDirectory] = getCommandOption(options_key[OptionKeyword::ReplicatedNodeDirectory]);
}