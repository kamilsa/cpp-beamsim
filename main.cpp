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

NS_LOG_COMPONENT_DEFINE("GossipApp");

using namespace ns3;

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
            // Print status (less frequently to reduce output)
            // if (m_receivedSignatures.size() % 100 == 0) {
            //     PrintWithTime("Peer " + std::to_string(m_peerId) +
            //                   " has received " + std::to_string(m_receivedSignatures.size()) +
            //                   " unique signatures");
            // }

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
    }

    void ReceiveSubnetAggregation(uint32_t subnetId) {
        if (m_receivedSubnets.insert(subnetId).second) {
            PrintWithTime("GlobalAggregator received aggregation from subnet " + std::to_string(subnetId));
            if (m_receivedSubnets.size() == m_nSubnets) {
                PrintWithTime(
                    "All subnet aggregations received (" + std::to_string(m_nSubnets) + "/" + std::to_string(m_nSubnets)
                    + "). Stopping simulation.");
                Simulator::Stop();
            }
        }
    }

private:
    uint32_t m_nSubnets;
    std::set<uint32_t> m_receivedSubnets;

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
               std::vector<Ptr<PeerApp> > *peerApps = nullptr) {
        m_subnetId = subnetId;
        m_nPeers = nPeers;
        m_globalAgg = globalAgg;
        m_peerApps = peerApps;
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
        m_globalAgg->ReceiveSubnetAggregation(m_subnetId);
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

    void PrintWithTime(const std::string &msg) {
        double timeInSeconds = Simulator::Now().GetSeconds();
        int milliseconds = static_cast<int>((timeInSeconds - std::floor(timeInSeconds)) * 1000);
        std::cout << std::fixed << std::setprecision(0) << std::floor(timeInSeconds) << "."
                << std::setfill('0') << std::setw(3) << milliseconds << " ms: " << msg << std::endl;
    }
};

int main(int argc, char *argv[]) {
    LogComponentEnable("GossipApp", LOG_LEVEL_INFO);

    // Record start time
    std::chrono::steady_clock::time_point realStartTime = std::chrono::steady_clock::now();

    const uint32_t nSubnets = 5;
    const uint32_t nPeersPerSubnet = 256;
    NodeContainer globalAggNode;
    globalAggNode.Create(1);
    Ptr<GlobalAggregatorApp> globalAggApp = CreateObject<GlobalAggregatorApp>();
    globalAggApp->Setup(nSubnets);
    globalAggNode.Get(0)->AddApplication(globalAggApp);
    std::vector<Ptr<SubnetAggregatorApp> > subnetAggApps(nSubnets);
    std::vector<NodeContainer> subnetNodes(nSubnets);
    std::vector<std::vector<Ptr<PeerApp> > > peerApps(nSubnets);
    // Create subnets and install apps
    for (uint32_t s = 0; s < nSubnets; ++s) {
        subnetNodes[s].Create(nPeersPerSubnet + 1); // +1 for aggregator
        Ptr<SubnetAggregatorApp> aggApp = CreateObject<SubnetAggregatorApp>();
        aggApp->Setup(s, nPeersPerSubnet, globalAggApp, &peerApps[s]);
        subnetNodes[s].Get(nPeersPerSubnet)->AddApplication(aggApp);
        subnetAggApps[s] = aggApp;
        peerApps[s].resize(nPeersPerSubnet);
        // Install PeerApps
        for (uint32_t p = 0; p < nPeersPerSubnet; ++p) {
            Ptr<PeerApp> peerApp = CreateObject<PeerApp>();
            peerApp->Setup(p, aggApp, &peerApps[s], nPeersPerSubnet);
            subnetNodes[s].Get(p)->AddApplication(peerApp);
            peerApps[s][p] = peerApp;
        }
    }

    Simulator::Stop(Seconds(10.0)); // Failsafe
    std::cout << "Starting ns-3 subnet aggregation simulation" << std::endl;

    // Record simulation start time
    double simStartTime = Simulator::Now().GetSeconds();

    Simulator::Run();

    // Get simulation end time
    double simEndTime = Simulator::Now().GetSeconds();
    // Get real end time
    std::chrono::steady_clock::time_point realEndTime = std::chrono::steady_clock::now();

    // Calculate time differences
    double virtualTimeSeconds = simEndTime - simStartTime;
    double realTimeSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(realEndTime - realStartTime).count() / 1000.0;

    // Print statistics
    std::cout << "\n--- SIMULATION STATISTICS ---" << std::endl;
    std::cout << "Total virtual time: " << std::fixed << std::setprecision(3) << virtualTimeSeconds << " seconds" << std::endl;
    std::cout << "Total real execution time: " << std::fixed << std::setprecision(3) << realTimeSeconds << " seconds" << std::endl;
    std::cout << "Ratio (virtual/real): " << std::fixed << std::setprecision(6) << (virtualTimeSeconds / realTimeSeconds) << std::endl;
    std::cout << "Total number of peers: " << nSubnets * nPeersPerSubnet << std::endl;
    std::cout << "----------------------------" << std::endl;

    Simulator::Destroy();
    return 0;
}
