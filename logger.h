#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <ns3/simulator.h>

namespace BeamSimLogger {

// Log levels
enum LogLevel {
    INFO,
    WARNING,
    ERROR,
    ROUND_CHANGE
};

// Centralized logging function
inline void Log(const std::string& message, LogLevel level = INFO) {
    double timeInSeconds = ns3::Simulator::Now().GetSeconds();
    int milliseconds = static_cast<int>((timeInSeconds - std::floor(timeInSeconds)) * 1000);

    // Start with the timestamp
    std::cout << std::fixed << std::setprecision(0) << std::floor(timeInSeconds)
              << "." << std::setfill('0') << std::setw(3) << milliseconds << " ms: ";

    // For round changes, add special formatting
    if (level == ROUND_CHANGE) {
        std::cout << "\n==========================================" << std::endl;
        std::cout << "==== " << message << " ====" << std::endl;
        std::cout << "==========================================" << std::endl;
    }
    // Regular log entry
    else {
        std::cout << message << std::endl;
    }
}

// Specialized log for round information
inline void LogRound(const std::string& message, int round, int totalRounds) {
    Log(message + " [Round " + std::to_string(round) + "/" + std::to_string(totalRounds) + "]");
}

// Specialized log for MPI process information
inline void LogMPI(const std::string& message, int process) {
    Log(message + " [Process " + std::to_string(process) + "]");
}

// Specialized log for combined round and process information
inline void LogRoundMPI(const std::string& message, int round, int totalRounds, int process) {
    Log(message + " [Round " + std::to_string(round) + "/" + std::to_string(totalRounds) +
        ", Process " + std::to_string(process) + "]");
}

// Log for round changes
inline void LogRoundChange(int previousRound, int newRound, int totalRounds) {
    Log("Round " + std::to_string(previousRound) + " complete. Starting round " +
        std::to_string(newRound) + " of " + std::to_string(totalRounds), ROUND_CHANGE);
}

} // namespace BeamSimLogger

#endif // LOGGER_H
