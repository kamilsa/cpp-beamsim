// src/main.cc
#include <ns3/core-module.h>
#include <ns3/network-module.h>
#include <ns3/internet-module.h>
#include <set>
#include <vector>
#include <algorithm>
#include <random>
#include <map>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <mpi.h>  // Include MPI header

NS_LOG_COMPONENT_DEFINE("GossipApp");

using namespace ns3;

// Structure to hold signature information for MPI communication
struct SignatureMessage {
    uint32_t originPeerId;
    uint32_t subnetId;

    SignatureMessage(uint32_t origin = 0, uint32_t subnet = 0)
        : originPeerId(origin), subnetId(subnet) {}
};

// --- PeerApp ---
class SubnetAggregatorApp;

class PeerApp : public Application {
public:
    void Setup(uint32_t peerId, Ptr<SubnetAggregatorApp> aggregator, std::vector<Ptr<PeerApp>>* peerApps = nullptr, uint32_t nPeers = 0) {
        m_peerId = peerId;
        m_aggregator = aggregator;
        m_peerApps = peerApps;
        m_nPeers = nPeers;
        m_fanOut = 6;  // Gossip fan-out factor
    }

    virtual void StartApplication() override {
        // Generate own signature as the originator
        Simulator::Schedule(Seconds(0.01), &PeerApp::GossipSignature, this, m_peerId);
    }

    void GossipSignature(uint32_t originPeerId) {
        // If we haven't seen this signature before, process and forward it
        if (m_receivedSignatures.insert(originPeerId).second) {
            // Select random peers to forward to (fan-out peers)
            if (m_peerApps && m_peerApps->size() > 0) {
                // Generate random targets for gossip
                std::vector<uint32_t> targets = SelectRandomPeers(m_fanOut);

                // Forward to selected peers
                for (uint32_t targetId: targets) {
                    // Avoid sending to self or to original source
                    if (targetId != m_peerId && targetId != originPeerId) {
                        // Add small delay for each gossip message
                        double delay = 0.001 + (0.001 * ((double) rand() / RAND_MAX)); // Small random delay
                        Simulator::Schedule(Seconds(delay), &PeerApp::ReceiveGossipSignature,
                                            (*m_peerApps)[targetId], originPeerId);
                    }
                }
            }
        }
    }

    void ReceiveGossipSignature(uint32_t originPeerId) {
        // Process the gossip message
        GossipSignature(originPeerId);
    }

    // This is kept for backward compatibility
    void SendSignature() {
        // Now handled by GossipSignature
    }

    // Method to check if peer has received a signature
    bool HasReceivedSignature(uint32_t originPeerId) const {
        return m_receivedSignatures.find(originPeerId) != m_receivedSignatures.end();
    }

    // Getter for peer ID
    uint32_t GetPeerId() const {
        return m_peerId;
    }

    // Get all received signatures
    const std::set<uint32_t> &GetReceivedSignatures() const {
        return m_receivedSignatures;
    }

private:
    uint32_t m_peerId;
    uint32_t m_nPeers;
    uint32_t m_fanOut;
    Ptr<SubnetAggregatorApp> m_aggregator;
    std::set<uint32_t> m_receivedSignatures; // Track signatures we've seen
    std::vector<Ptr<PeerApp> > *m_peerApps; // Pointer to all peer apps in the subnet

    void PrintWithTime(const std::string &msg) {
        double timeInSeconds = Simulator::Now().GetSeconds();
        int milliseconds = static_cast<int>((timeInSeconds - std::floor(timeInSeconds)) * 1000);
        std::cout << std::fixed << std::setprecision(0) << std::floor(timeInSeconds) << "."
                << std::setfill('0') << std::setw(3) << milliseconds << " ms: " << msg << std::endl;
    }

    std::vector<uint32_t> SelectRandomPeers(uint32_t count) {
        std::vector<uint32_t> result;
        std::vector<uint32_t> candidates;

        // Create list of candidate peers (all except self)
        for (uint32_t i = 0; i < m_nPeers; i++) {
            if (i != m_peerId) {
                candidates.push_back(i);
            }
        }

        // If we have fewer candidates than fan-out, use all of them
        uint32_t selectCount = std::min(count, static_cast<uint32_t>(candidates.size()));

        // Shuffle and select the first 'selectCount' peers
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(candidates.begin(), candidates.end(), g);

        for (uint32_t i = 0; i < selectCount; i++) {
            result.push_back(candidates[i]);
        }

        return result;
    }
};

// --- GlobalAggregatorApp ---
class GlobalAggregatorApp : public Application {
public:
    void Setup(uint32_t nSubnets) {
        m_nSubnets = nSubnets;
        m_completionTime = 0.0;
    }

    void ReceiveSubnetAggregation(uint32_t subnetId) {
        if (m_receivedSubnets.insert(subnetId).second) {
            PrintWithTime("GlobalAggregator received aggregation from subnet " + std::to_string(subnetId));
            if (m_receivedSubnets.size() == m_nSubnets) {
                m_completionTime = Simulator::Now().GetSeconds(); // Record when work actually completes
                PrintWithTime(
                    "All subnet aggregations received (" + std::to_string(m_nSubnets) + "/" + std::to_string(m_nSubnets)
                    + "). Stopping simulation.");
                Simulator::Stop();
            }
        }
    }

    double GetCompletionTime() const {
        return m_completionTime;
    }

private:
    uint32_t m_nSubnets;
    std::set<uint32_t> m_receivedSubnets;
    double m_completionTime; // Store the actual completion time

    void PrintWithTime(const std::string &msg) {
        double timeInSeconds = Simulator::Now().GetSeconds();
        int milliseconds = static_cast<int>((timeInSeconds - std::floor(timeInSeconds)) * 1000);
        std::cout << std::fixed << std::setprecision(0) << std::floor(timeInSeconds) << "."
                << std::setfill('0') << std::setw(3) << milliseconds << " ms: " << msg << std::endl;
    }
};

// --- SubnetAggregatorApp ---
class SubnetAggregatorApp : public Application {
public:
    void Setup(uint32_t subnetId, uint32_t nPeers, Ptr<GlobalAggregatorApp> globalAgg,
               std::vector<Ptr<PeerApp> > *peerApps = nullptr, int mpiRank = 0) {
        m_subnetId = subnetId;
        m_nPeers = nPeers;
        m_globalAgg = globalAgg;
        m_peerApps = peerApps;
        m_mpiRank = mpiRank;
        m_mpiMaster = (mpiRank == 0);
    }

    void ReceiveSignature(uint32_t peerId) {
        if (m_receivedPeers.insert(peerId).second) {
            if (m_receivedPeers.size() % 10 == 0) {
                PrintWithTime("SA from subnet " + std::to_string(m_subnetId) +
                              " received " + std::to_string(m_receivedPeers.size()) + " signatures from subnet " +
                              std::to_string(m_subnetId));
            }
            if (m_receivedPeers.size() == m_threshold) {
                PrintWithTime("SA from subnet " + std::to_string(m_subnetId) +
                              " reached threshold of " + std::to_string(m_threshold) +
                              " signatures. Preparing to send aggregation.");
                Simulator::Schedule(Seconds(0.01), &SubnetAggregatorApp::SendAggregation, this);
            }
        }
    }

    virtual void StartApplication() override {
        m_threshold = m_nPeers * 2 / 3 + 1;

        // Schedule periodic monitoring of signatures spread through gossip
        Simulator::Schedule(Seconds(0.1), &SubnetAggregatorApp::MonitorSignatures, this);
    }

    void SendAggregation() {
        PrintWithTime("SA from subnet " + std::to_string(m_subnetId) +
                      " sending aggregation to global aggregator");

        if (m_mpiMaster) {
            // Only send to global aggregator if we're on the master process
            m_globalAgg->ReceiveSubnetAggregation(m_subnetId);
        } else {
            // Send to master process via MPI
            int destination = 0; // Master rank
            int tag = 1000 + m_subnetId; // Unique tag for each subnet
            uint32_t data = m_subnetId;
            MPI_Send(&data, 1, MPI_UINT32_T, destination, tag, MPI_COMM_WORLD);

            PrintWithTime("SA from subnet " + std::to_string(m_subnetId) +
                          " sent aggregation to master process via MPI");
        }
    }

    void MonitorSignatures() {
        // Only proceed if we have peers to monitor
        if (m_peerApps && m_peerApps->size() > 0) {
            // Check each peer's signatures and update our record
            for (uint32_t i = 0; i < m_nPeers; i++) {
                Ptr<PeerApp> peer = (*m_peerApps)[i];
                const std::set<uint32_t> &peerSigs = peer->GetReceivedSignatures();

                // Add all signatures from this peer to our record
                for (uint32_t sigId: peerSigs) {
                    ReceiveSignature(sigId);
                }
            }
        }

        // Schedule next monitoring round
        Simulator::Schedule(Seconds(0.1), &SubnetAggregatorApp::MonitorSignatures, this);
    }

private:
    uint32_t m_subnetId;
    uint32_t m_nPeers;
    uint32_t m_threshold;
    std::set<uint32_t> m_receivedPeers;
    Ptr<GlobalAggregatorApp> m_globalAgg;
    std::vector<Ptr<PeerApp> > *m_peerApps; // Pointer to all peer apps in the subnet
    int m_mpiRank;
    bool m_mpiMaster;

    void PrintWithTime(const std::string &msg) {
        double timeInSeconds = Simulator::Now().GetSeconds();
        int milliseconds = static_cast<int>((timeInSeconds - std::floor(timeInSeconds)) * 1000);
        std::cout << std::fixed << std::setprecision(0) << std::floor(timeInSeconds) << "."
                << std::setfill('0') << std::setw(3) << milliseconds << " ms: " << msg << std::endl;
    }
};

// Receive messages from worker processes (in master process)
void CheckForMPIMessages(Ptr<GlobalAggregatorApp> globalAggApp, int world_size) {
    MPI_Status status;
    int flag = 0;
    uint32_t data;

    // Check for messages from any source
    for (int source = 1; source < world_size; source++) {
        MPI_Iprobe(source, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

        if (flag) {
            // Receive the message
            MPI_Recv(&data, 1, MPI_UINT32_T, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Process the subnet aggregation
            globalAggApp->ReceiveSubnetAggregation(data);

            std::cout << Simulator::Now().GetSeconds() << " ms: Master received subnet aggregation via MPI from subnet "
                      << data << " (process " << status.MPI_SOURCE << ")" << std::endl;
        }
    }

    // Schedule next check
    Simulator::Schedule(Seconds(0.01), &CheckForMPIMessages, globalAggApp, world_size);
}

int main(int argc, char *argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get MPI process information
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Create MPI datatype for SignatureMessage
    MPI_Datatype mpi_signature_type;
    int blocklengths[2] = {1, 1};
    MPI_Aint displacements[2] = {0, sizeof(uint32_t)};
    MPI_Datatype types[2] = {MPI_UINT32_T, MPI_UINT32_T};
    MPI_Type_create_struct(2, blocklengths, displacements, types, &mpi_signature_type);
    MPI_Type_commit(&mpi_signature_type);

    LogComponentEnable("GossipApp", LOG_LEVEL_INFO);

    // Record start time
    std::chrono::steady_clock::time_point realStartTime = std::chrono::steady_clock::now();

    // Configuration
    const uint32_t nSubnets = 100; // Total number of subnets across all processes
    const uint32_t nPeersPerSubnet = 1024;

    // Calculate which subnets to process on this rank
    uint32_t subnets_per_rank = nSubnets / world_size;
    uint32_t extra_subnets = nSubnets % world_size;
    uint32_t my_subnet_start = world_rank * subnets_per_rank + std::min(world_rank, static_cast<int>(extra_subnets));
    uint32_t my_subnet_count = subnets_per_rank + (world_rank < static_cast<int>(extra_subnets) ? 1 : 0);
    uint32_t my_subnet_end = my_subnet_start + my_subnet_count;

    if (world_rank == 0) {
        std::cout << "Running simulation with " << world_size << " MPI processes" << std::endl;
        std::cout << "Total subnets: " << nSubnets << ", peers per subnet: " << nPeersPerSubnet << std::endl;
    }

    std::cout << "Process " << world_rank << " handling subnets " << my_subnet_start
              << " to " << (my_subnet_end-1) << " (" << my_subnet_count << " subnets)" << std::endl;

    // Global aggregator is needed only on the master process (rank 0)
    Ptr<GlobalAggregatorApp> globalAggApp;
    NodeContainer globalAggNode;

    if (world_rank == 0) {
        globalAggNode.Create(1);
        globalAggApp = CreateObject<GlobalAggregatorApp>();
        globalAggApp->Setup(nSubnets);
        globalAggNode.Get(0)->AddApplication(globalAggApp);

        // Schedule regular checks for MPI messages
        Simulator::Schedule(Seconds(0.1), &CheckForMPIMessages, globalAggApp, world_size);
    }

    // Each process creates only the subnets it's responsible for
    std::vector<Ptr<SubnetAggregatorApp>> subnetAggApps(my_subnet_count);
    std::vector<NodeContainer> subnetNodes(my_subnet_count);
    std::vector<std::vector<Ptr<PeerApp>>> peerApps(my_subnet_count);

    // Create only subnets assigned to this process
    for (uint32_t i = 0; i < my_subnet_count; i++) {
        uint32_t s = my_subnet_start + i;  // Global subnet ID
        subnetNodes[i].Create(nPeersPerSubnet + 1);  // +1 for aggregator

        Ptr<SubnetAggregatorApp> aggApp = CreateObject<SubnetAggregatorApp>();
        aggApp->Setup(s, nPeersPerSubnet, globalAggApp, &peerApps[i], world_rank);
        subnetNodes[i].Get(nPeersPerSubnet)->AddApplication(aggApp);
        subnetAggApps[i] = aggApp;

        peerApps[i].resize(nPeersPerSubnet);

        // Install PeerApps
        for (uint32_t p = 0; p < nPeersPerSubnet; p++) {
            Ptr<PeerApp> peerApp = CreateObject<PeerApp>();
            peerApp->Setup(p, aggApp, &peerApps[i], nPeersPerSubnet);
            subnetNodes[i].Get(p)->AddApplication(peerApp);
            peerApps[i][p] = peerApp;
        }
    }

    Simulator::Stop(Seconds(10.0));  // Failsafe

    if (world_rank == 0) {
        std::cout << "Starting ns-3 subnet aggregation simulation" << std::endl;
    }

    // Record simulation start time
    double simStartTime = Simulator::Now().GetSeconds();

    // Create a barrier to ensure all processes start simulation at roughly the same time
    MPI_Barrier(MPI_COMM_WORLD);

    // Run the simulation
    Simulator::Run();

    // Synchronize all processes before collecting results
    MPI_Barrier(MPI_COMM_WORLD);

    // Get simulation end time and actual completion time
    double simEndTime = Simulator::Now().GetSeconds();

    // Get the actual completion time from the global aggregator on master process
    double actualCompletionTime = 0.0;
    if (world_rank == 0 && globalAggApp) {
        actualCompletionTime = globalAggApp->GetCompletionTime();
        // If completion time is 0, it means simulation terminated at failsafe
        if (actualCompletionTime < 0.001) {
            actualCompletionTime = simEndTime;  // Use the simulator end time instead
        }
    }

    // Broadcast the actual completion time to all processes
    MPI_Bcast(&actualCompletionTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Get real end time
    std::chrono::steady_clock::time_point realEndTime = std::chrono::steady_clock::now();

    // Calculate time differences - use actual completion time for virtual time
    double virtualTimeSeconds = actualCompletionTime - simStartTime;
    double realTimeSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(realEndTime - realStartTime).count() / 1000.0;

    // Collect timing statistics from all processes
    double max_virtualTimeSeconds;
    double max_realTimeSeconds;
    MPI_Reduce(&virtualTimeSeconds, &max_virtualTimeSeconds, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&realTimeSeconds, &max_realTimeSeconds, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Print statistics on the master process
    if (world_rank == 0) {
        std::cout << "\n--- SIMULATION STATISTICS ---" << std::endl;
        std::cout << "MPI processes: " << world_size << std::endl;
        std::cout << "Total virtual time: " << std::fixed << std::setprecision(3) << max_virtualTimeSeconds << " seconds" << std::endl;
        std::cout << "Total real execution time: " << std::fixed << std::setprecision(3) << max_realTimeSeconds << " seconds" << std::endl;
        std::cout << "Ratio (virtual/real): " << std::fixed << std::setprecision(6) << (max_virtualTimeSeconds / max_realTimeSeconds) << std::endl;
        std::cout << "Total number of peers: " << nSubnets * nPeersPerSubnet << std::endl;
        std::cout << "----------------------------" << std::endl;
    }

    Simulator::Destroy();

    // Clean up MPI datatype
    MPI_Type_free(&mpi_signature_type);

    MPI_Finalize();  // Finalize MPI
    return 0;
}
