// src/main.cc
#include <ns3/core-module.h>
#include <ns3/network-module.h>
#include <ns3/internet-module.h>
#include <set>
#include <vector>
#include <algorithm>
#include <random>
#include <map>

NS_LOG_COMPONENT_DEFINE("GossipApp");

using namespace ns3;

// --- PeerApp ---
class SubnetAggregatorApp;

class PeerApp : public Application {
public:
    void Setup(uint32_t peerId, Ptr<SubnetAggregatorApp> aggregator) {
        m_peerId = peerId;
        m_aggregator = aggregator;
    }

    virtual void StartApplication() override {
        Simulator::Schedule(Seconds(0.01), &PeerApp::SendSignature, this);
    }

    void SendSignature();

private:
    uint32_t m_peerId;
    Ptr<SubnetAggregatorApp> m_aggregator;
    void PrintWithTime(const std::string& msg) {
        std::cout << Simulator::Now().GetSeconds() << ": " << msg << std::endl;
    }
};

// --- SubnetAggregatorApp ---
class GlobalAggregatorApp;

class SubnetAggregatorApp : public Application {
public:
    void Setup(uint32_t subnetId, uint32_t nPeers, Ptr<GlobalAggregatorApp> globalAgg) {
        m_subnetId = subnetId;
        m_nPeers = nPeers;
        m_globalAgg = globalAgg;
    }

    void ReceiveSignature(uint32_t peerId) {
        if (m_receivedPeers.insert(peerId).second) {
            if (m_receivedPeers.size() % 10 == 0) {
                std::cout << Simulator::Now().GetSeconds() << ": SubnetAggregator from subnet " << m_subnetId <<
                        " received " << m_receivedPeers.size() << " signatures from subnet " <<
                        m_subnetId << std::endl;
            }
            if (m_receivedPeers.size() == m_threshold) {
                Simulator::Schedule(Seconds(0.01), &SubnetAggregatorApp::SendAggregation, this);
            }
        }
    }

    virtual void StartApplication() override {
        m_threshold = m_nPeers * 2 / 3 + 1;
    }

    void SendAggregation();

private:
    uint32_t m_subnetId;
    uint32_t m_nPeers;
    uint32_t m_threshold;
    std::set<uint32_t> m_receivedPeers;
    Ptr<GlobalAggregatorApp> m_globalAgg;
};

// --- GlobalAggregatorApp ---
class GlobalAggregatorApp : public Application {
public:
    void Setup(uint32_t nSubnets) {
        m_nSubnets = nSubnets;
    }

    void ReceiveSubnetAggregation(uint32_t subnetId) {
        if (m_receivedSubnets.insert(subnetId).second) {
            std::cout << Simulator::Now().GetSeconds() << ": GlobalAggregator received aggregation from subnet " << subnetId << std::endl;
            if (m_receivedSubnets.size() == m_nSubnets) {
                std::cout << Simulator::Now().GetSeconds() << ": All subnet aggregations received. Stopping simulation." << std::endl;
                Simulator::Stop();
            }
        }
    }

private:
    uint32_t m_nSubnets;
    std::set<uint32_t> m_receivedSubnets;
};

// --- Implementation of PeerApp::SendSignature and SubnetAggregatorApp::SendAggregation ---
void PeerApp::SendSignature() {
    m_aggregator->ReceiveSignature(m_peerId);
}

void SubnetAggregatorApp::SendAggregation() {
    m_globalAgg->ReceiveSubnetAggregation(m_subnetId);
}

int main(int argc, char *argv[]) {
    LogComponentEnable("GossipApp", LOG_LEVEL_INFO);
    NS_LOG_INFO("Starting ns-3 subnet aggregation simulation");
    const uint32_t nSubnets = 10;
    const uint32_t nPeersPerSubnet = 1024;
    NodeContainer globalAggNode;
    globalAggNode.Create(1);
    Ptr<GlobalAggregatorApp> globalAggApp = CreateObject<GlobalAggregatorApp>();
    globalAggApp->Setup(nSubnets);
    globalAggNode.Get(0)->AddApplication(globalAggApp);
    std::vector<Ptr<SubnetAggregatorApp> > subnetAggApps(nSubnets);
    std::vector<NodeContainer> subnetNodes(nSubnets);
    // Create subnets and install apps
    for (uint32_t s = 0; s < nSubnets; ++s) {
        subnetNodes[s].Create(nPeersPerSubnet + 1); // +1 for aggregator
        Ptr<SubnetAggregatorApp> aggApp = CreateObject<SubnetAggregatorApp>();
        aggApp->Setup(s, nPeersPerSubnet, globalAggApp);
        subnetNodes[s].Get(nPeersPerSubnet)->AddApplication(aggApp);
        subnetAggApps[s] = aggApp;
        // Install PeerApps
        for (uint32_t p = 0; p < nPeersPerSubnet; ++p) {
            Ptr<PeerApp> peerApp = CreateObject<PeerApp>();
            peerApp->Setup(p, aggApp);
            subnetNodes[s].Get(p)->AddApplication(peerApp);
        }
    }
    Simulator::Stop(Seconds(10.0)); // Failsafe
    Simulator::Run();
    Simulator::Destroy();
    return 0;
}
