#ifndef COMMANDLIST_H
#define COMMANDLIST_H

#include "utility/CommandParser.h"
#include <map>
#include <iostream>

enum OptionKeyword {
    QueryGraphFile = 0,                     // -q, The query graph file path
    PatitionedDataGraphDirectory = 1,       // -d, The partitioned data graph directory path
    ThreadCount = 2,                        // -n, The number of thread used in thread pool, default is 4
    k = 3,                                  // -k, The top-k value
    DataGraphFile = 4,                      // -g, The data graph file path
    PartitionCount = 5,                     // -p, Number of partitions
    DistanceFile = 6,                       // -m, The distance matrix of PDG
    GroundTruthFile = 7,                    // -t, The ground truth file of partition prediction stage
    ReplicatedNodeDirectory = 8,            // -r, The directory of files recording nodes in all partitions which are replicated during parititoning
};

class MatchingCommand : public CommandParser{
private:
    std::map<OptionKeyword, std::string> options_key;
    std::map<OptionKeyword, std::string> options_value;

private:
    void processOptions();

public:
    MatchingCommand(int argc, char **argv);

    std::string getQueryGraphFilePath() {
        return options_value[OptionKeyword::QueryGraphFile];
    }

    std::string getPartitionedDataGraphDirectoryPath() {
        return options_value[OptionKeyword::PatitionedDataGraphDirectory];
    }

    std::string getThreadCount() {
        return options_value[OptionKeyword::ThreadCount] == "" ? "1" : options_value[OptionKeyword::ThreadCount];
    }

    std::string getK() {
        return options_value[OptionKeyword::k] == "" ? "None" : options_value[OptionKeyword::k];
    }

    std::string getDataGraphFilePath() {
        return options_value[OptionKeyword::DataGraphFile];
    }

    std::string getPartitionNum() {
        return options_value[OptionKeyword::PartitionCount] == "" ? "None" : options_value[OptionKeyword::PartitionCount];
    }

    std::string getDistanceFilePath() {
        return options_value[OptionKeyword::DistanceFile];
    }

    std::string getGroundTruthFilePath() {
        return options_value[OptionKeyword::GroundTruthFile];
    }

    std::string getReplicatedNodeDirectory() {
        return options_value[OptionKeyword::ReplicatedNodeDirectory];
    }
};

#endif //COMMANDLIST_H