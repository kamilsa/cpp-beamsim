// src/main.cc
#include <ns3/core-module.h>
#include <ns3/network-module.h>
#include <ns3/internet-module.h>
#include <ns3/point-to-point-module.h>
#include <ns3/applications-module.h>
#include <ns3/command-line.h>
#include <set>
#include <vector>
#include <algorithm>
#include <random>
#include <map>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <mpi.h>  // Include MPI header
#include "logger.h"  // Include our custom logger

NS_LOG_COMPONENT_DEFINE("GossipApp");

using namespace ns3;
using namespace BeamSimLogger;  // Using our logger namespace

// Forward declarations
class PeerApp;

// Define a custom packet type for signatures
class SignatureHeader : public Header
{
public:
  SignatureHeader();
  virtual ~SignatureHeader();

  void SetOriginPeerId(uint32_t originPeerId);
  uint32_t GetOriginPeerId(void) const;

  static TypeId GetTypeId(void);
  virtual TypeId GetInstanceTypeId(void) const;
  virtual void Print(std::ostream &os) const;
  virtual void Serialize(Buffer::Iterator start) const;
  virtual uint32_t Deserialize(Buffer::Iterator start);
  virtual uint32_t GetSerializedSize(void) const;

private:
  uint32_t m_originPeerId;
  // Add padding to simulate real signature size (1.5KB)
  uint8_t m_padding[1528];
};

SignatureHeader::SignatureHeader()
  : m_originPeerId(0)
{
  memset(m_padding, 0, sizeof(m_padding));
}

SignatureHeader::~SignatureHeader()
{
}

TypeId
SignatureHeader::GetTypeId(void)
{
  static TypeId tid = TypeId("ns3::SignatureHeader")
    .SetParent<Header>()
    .AddConstructor<SignatureHeader>();
  return tid;
}

TypeId
SignatureHeader::GetInstanceTypeId(void) const
{
  return GetTypeId();
}

void
SignatureHeader::Print(std::ostream &os) const
{
  os << "OriginPeerId=" << m_originPeerId;
}

uint32_t
SignatureHeader::GetSerializedSize(void) const
{
  return sizeof(m_originPeerId) + sizeof(m_padding); // 4 bytes + padding bytes
}

void
SignatureHeader::Serialize(Buffer::Iterator start) const
{
  start.WriteHtonU32(m_originPeerId);
  start.Write(m_padding, sizeof(m_padding));
}

uint32_t
SignatureHeader::Deserialize(Buffer::Iterator start)
{
  m_originPeerId = start.ReadNtohU32();
  start.Read(m_padding, sizeof(m_padding));
  return GetSerializedSize();
}

void
SignatureHeader::SetOriginPeerId(uint32_t originPeerId)
{
  m_originPeerId = originPeerId;
}

uint32_t
SignatureHeader::GetOriginPeerId(void) const
{
  return m_originPeerId;
}

// Signature Server (receiver)
class SignatureServer : public Application
{
public:
  static TypeId GetTypeId(void);
  SignatureServer();
  virtual ~SignatureServer();

  void SetPeerApp(Ptr<PeerApp> peerApp);

protected:
  virtual void DoDispose(void);

private:
  virtual void StartApplication(void);
  virtual void StopApplication(void);

  void HandleRead(Ptr<Socket> socket);

  Ptr<Socket> m_socket;
  uint16_t m_port;
  Ptr<PeerApp> m_peerApp; // Reference to the parent PeerApp
};

// Signature Client (sender)
class SignatureClient : public Application
{
public:
  static TypeId GetTypeId(void);
  SignatureClient();
  virtual ~SignatureClient();

  void SetAttribute(std::string name, const AttributeValue &value);
  void SendSignature(uint32_t originPeerId);
  void Setup(Ipv4Address address, uint16_t port);

protected:
  virtual void DoDispose(void);

private:
  virtual void StartApplication(void);
  virtual void StopApplication(void);

  void ConnectionSucceeded(Ptr<Socket> socket);
  void ConnectionFailed(Ptr<Socket> socket);

  Ptr<Socket> m_socket;
  Address m_peerAddress;
  uint16_t m_peerPort;
  bool m_connected;
  uint32_t m_lastOriginPeerId; // For retrying if needed
};

// Structure to hold signature information for MPI communication
struct SignatureMessage {
    uint32_t originPeerId;
    uint32_t subnetId;
    // Add padding to make the signature size 1.5KB (1536 bytes)
    // 1536 - 8 (two uint32_t fields) = 1528 bytes of padding
    char signatureData[1528];

    SignatureMessage(uint32_t origin = 0, uint32_t subnet = 0)
        : originPeerId(origin), subnetId(subnet) {
        // Initialize with zeros - no need for random data in simulation
        std::memset(signatureData, 0, sizeof(signatureData));
    }
};

// --- PeerApp ---
class SubnetAggregatorApp;

class PeerApp : public Application {
public:
    // Add static counter for total messages sent
    static uint64_t s_totalMessagesSent;

    void Setup(uint32_t peerId, Ptr<SubnetAggregatorApp> aggregator, std::vector<Ptr<PeerApp>>* peerApps = nullptr, uint32_t nPeers = 0) {
        m_peerId = peerId;
        m_aggregator = aggregator;
        m_peerApps = peerApps;
        m_nPeers = nPeers;
        m_fanOut = 6;  // Gossip fan-out factor
        m_currentRound = 1;

        // For P2P communication
        m_signatureClients.clear();
        m_signatureServer = nullptr;
    }

    // Set up P2P connections with other peers
    void SetupP2PConnections(const std::vector<Ipv4Address>& peerAddresses) {
        m_peerAddresses = peerAddresses;
    }

    // Add connection to a specific peer
    void AddPeerConnection(uint32_t peerId, Ipv4Address address) {
        if (peerId == m_peerId) return;  // Don't connect to self

        if (m_connectedPeers.find(peerId) == m_connectedPeers.end()) {
            m_connectedPeers.insert(peerId);

            // Create a client for this connection
            Ptr<SignatureClient> client = CreateObject<SignatureClient>();
            client->Setup(address, 9000);  // Port 9000 for all signature servers
            GetNode()->AddApplication(client);
            m_signatureClients[peerId] = client;
        }
    }

    // Set the local server for receiving signatures
    void SetSignatureServer(Ptr<SignatureServer> server) {
        m_signatureServer = server;
        m_signatureServer->SetPeerApp(this);
    }

    virtual void StartApplication() override {
        // Start the signature server if it's set
        if (m_signatureServer) {
            m_signatureServer->SetStartTime(Seconds(0.0));
        }

        // Start all clients
        for (auto& client : m_signatureClients) {
            client.second->SetStartTime(Seconds(0.0));
        }

        // Generate own signature as the originator
        Simulator::Schedule(Seconds(0.01), &PeerApp::GossipSignature, this, m_peerId);
    }

    void GossipSignature(uint32_t originPeerId) {
        // If we haven't seen this signature before in this round, process and forward it
        if (m_receivedSignatures.insert(originPeerId).second) {
            // Select random peers to forward to (fan-out peers)
            if (m_peerApps && m_peerApps->size() > 0) {
                // Generate random targets for gossip - only select peers we have connections with
                std::vector<uint32_t> targets = SelectRandomConnectedPeers(m_fanOut);

                // Forward to selected peers
                for (uint32_t targetId : targets) {
                    // Avoid sending to self or to original source
                    if (targetId != m_peerId && targetId != originPeerId) {
                        // Send via P2P if we have a connection
                        if (m_signatureClients.find(targetId) != m_signatureClients.end()) {
                            // Add small delay for each gossip message
                            double delay = 0.001 + (0.001 * ((double) rand() / RAND_MAX));
                            Simulator::Schedule(Seconds(delay), &PeerApp::SendSignatureToClient,
                                              this, targetId, originPeerId);
                        }
                    }
                }
            }
        }
    }

    // Helper to send a signature to a specific peer via P2P
    void SendSignatureToClient(uint32_t targetPeerId, uint32_t originPeerId) {
        auto it = m_signatureClients.find(targetPeerId);
        if (it != m_signatureClients.end()) {
            it->second->SendSignature(originPeerId);
            // Message counting is now handled in SignatureClient::SendSignature
        }
    }

    void ReceiveGossipSignature(const SignatureMessage& sig) {
        // Process the gossip message - extract the origin peer ID
        GossipSignature(sig.originPeerId);
    }

    // Reset peer state for a new round
    void StartNewRound(uint32_t roundNumber) {
        m_currentRound = roundNumber;
        m_receivedSignatures.clear();

        // Start by sending own signature in the new round
        Simulator::Schedule(Seconds(0.01), &PeerApp::GossipSignature, this, m_peerId);
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

    // Static method to get total messages sent
    static uint64_t GetTotalMessagesSent() {
        return s_totalMessagesSent;
    }

private:
    uint32_t m_peerId;
    uint32_t m_nPeers;
    uint32_t m_fanOut;
    uint32_t m_currentRound;
    Ptr<SubnetAggregatorApp> m_aggregator;
    std::set<uint32_t> m_receivedSignatures; // Track signatures we've seen
    std::vector<Ptr<PeerApp> > *m_peerApps; // Pointer to all peer apps in the subnet

    // P2P communication
    std::vector<Ipv4Address> m_peerAddresses;
    std::map<uint32_t, Ptr<SignatureClient>> m_signatureClients;
    Ptr<SignatureServer> m_signatureServer;
    std::set<uint32_t> m_connectedPeers; // Peers we have connections to

    void PrintWithTime(const std::string &msg) {
        double timeInSeconds = Simulator::Now().GetSeconds();
        int milliseconds = static_cast<int>((timeInSeconds - std::floor(timeInSeconds)) * 1000);
        std::cout << std::fixed << std::setprecision(0) << std::floor(timeInSeconds) << "."
                << std::setfill('0') << std::setw(3) << milliseconds << " ms: " << msg << std::endl;
    }

    std::vector<uint32_t> SelectRandomPeers(uint32_t count) {
        // Original version - selects from all peers
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

    std::vector<uint32_t> SelectRandomConnectedPeers(uint32_t count) {
        // New version - only selects from peers we have connections to
        std::vector<uint32_t> result;
        std::vector<uint32_t> candidates;

        // Create list of connected peer IDs
        for (uint32_t peerId : m_connectedPeers) {
            if (peerId != m_peerId) {
                candidates.push_back(peerId);
            }
        }

        // If we have fewer candidates than fan-out, use all of them
        uint32_t selectCount = std::min(count, static_cast<uint32_t>(candidates.size()));

        if (selectCount > 0) {
            // Shuffle and select peers
            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(candidates.begin(), candidates.end(), g);

            for (uint32_t i = 0; i < selectCount; i++) {
                result.push_back(candidates[i]);
            }
        }

        return result;
    }
};

// Initialize static member
uint64_t PeerApp::s_totalMessagesSent = 0;

// --- GlobalAggregatorApp ---
class GlobalAggregatorApp : public Application {
public:
    void Setup(uint32_t nSubnets, uint32_t nRounds = 5) {
        m_nSubnets = nSubnets;
        m_completionTime = 0.0;
        m_currentRound = 1;
        m_totalRounds = nRounds;
    }

    void ReceiveSubnetAggregation(uint32_t subnetId) {
        if (m_receivedSubnets.insert(subnetId).second) {
            // Use the improved logging with consistent round information
            LogRound("GlobalAggregator received aggregation from subnet " + std::to_string(subnetId),
                     m_currentRound, m_totalRounds);

            if (m_receivedSubnets.size() == m_nSubnets) {
                if (m_currentRound < m_totalRounds) {
                    // Store previous round for logging
                    int prevRound = m_currentRound;

                    // Start next round
                    m_currentRound++;

                    // Use special formatting for round transitions with visual separator
                    LogRoundChange(prevRound, m_currentRound, m_totalRounds);

                    // Clear received subnets for the next round
                    m_receivedSubnets.clear();

                    // Notify all subnets to start next round
                    Simulator::Schedule(Seconds(0.01), &GlobalAggregatorApp::StartNextRound, this);
                } else {
                    // All rounds complete
                    m_completionTime = Simulator::Now().GetSeconds(); // Record when work actually completes
                    Log("All subnet aggregations received for all " + std::to_string(m_totalRounds) +
                        " rounds. Stopping simulation.");
                    Simulator::Stop();
                }
            }
        }
    }

    double GetCompletionTime() const {
        return m_completionTime;
    }

    uint32_t GetCurrentRound() const {
        return m_currentRound;
    }

    uint32_t GetTotalRounds() const {
        return m_totalRounds;
    }

    // Signal to start a new round
    void StartNextRound() {
        LogRound("GlobalAggregator signaling start of round " + std::to_string(m_currentRound),
                m_currentRound, m_totalRounds);
        // This event will be observed by the master MPI process which will broadcast to workers
    }

private:
    uint32_t m_nSubnets;
    uint32_t m_currentRound;
    uint32_t m_totalRounds;
    std::set<uint32_t> m_receivedSubnets;
    double m_completionTime; // Store the actual completion time
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
        m_currentRound = 1;
        m_totalRounds = globalAgg ? globalAgg->GetTotalRounds() : 3; // Default to 3 rounds if no global aggregator
    }

    void ReceiveSignature(uint32_t peerId) {
        if (m_receivedPeers.insert(peerId).second) {
            // if (m_receivedPeers.size() % 10 == 0) {
            //     LogRoundMPI("SA from subnet " + std::to_string(m_subnetId) +
            //                 " received " + std::to_string(m_receivedPeers.size()) + " signatures from subnet " +
            //                 std::to_string(m_subnetId),
            //                 m_currentRound, m_totalRounds, m_mpiRank);
            // }
            if (m_receivedPeers.size() == m_threshold) {
                LogRoundMPI("SA from subnet " + std::to_string(m_subnetId) +
                           " reached threshold of " + std::to_string(m_threshold) +
                           " signatures. Preparing to send aggregation",
                           m_currentRound, m_totalRounds, m_mpiRank);

                Simulator::Schedule(Seconds(0.01), &SubnetAggregatorApp::SendAggregation, this);
            }
        }
    }

    virtual void StartApplication() override {
        m_threshold = m_nPeers * 2 / 3 + 1;

        // Schedule periodic monitoring of signatures spread through gossip
        Simulator::Schedule(Seconds(0.1), &SubnetAggregatorApp::MonitorSignatures, this);

        // If we're on the master process, check for round changes
        if (m_mpiMaster) {
            Simulator::Schedule(Seconds(0.1), &SubnetAggregatorApp::CheckForNextRound, this);
        }
    }

    void SendAggregation() {
        LogRoundMPI("SA from subnet " + std::to_string(m_subnetId) +
                   " sending aggregation to global aggregator",
                   m_currentRound, m_totalRounds, m_mpiRank);

        if (m_mpiMaster) {
            // Only send to global aggregator if we're on the master process
            m_globalAgg->ReceiveSubnetAggregation(m_subnetId);
        } else {
            // Send to master process via MPI
            int destination = 0; // Master rank
            int tag = 1000 + m_subnetId; // Unique tag for each subnet
            uint32_t data = m_subnetId;
            MPI_Send(&data, 1, MPI_UINT32_T, destination, tag, MPI_COMM_WORLD);

            LogRoundMPI("SA from subnet " + std::to_string(m_subnetId) +
                       " sent aggregation to master process via MPI",
                       m_currentRound, m_totalRounds, m_mpiRank);
        }
    }

    void StartNextRound(uint32_t roundNumber) {
        // Reset for new round
        m_currentRound = roundNumber;
        m_receivedPeers.clear();

        // Notify peers to start a new round of gossip - but only if previous round is complete
        if (m_peerApps && m_peerApps->size() > 0) {
            LogRoundMPI("SA from subnet " + std::to_string(m_subnetId) +
                       " starting round " + std::to_string(m_currentRound),
                       m_currentRound, m_totalRounds, m_mpiRank);

            // Reset all peers for the new round first
            for (uint32_t i = 0; i < m_nPeers; i++) {
                Ptr<PeerApp> peer = (*m_peerApps)[i];
                peer->StartNewRound(m_currentRound);
            }

            // Initiate new round of gossip after a brief delay to ensure clean separation between rounds
            Simulator::Schedule(Seconds(0.05), &SubnetAggregatorApp::InitiateNewRound, this);
        }
    }

    void InitiateNewRound() {
        // Start the new round with each peer sending its own signature
        if (m_peerApps && m_peerApps->size() > 0) {
            for (uint32_t i = 0; i < m_nPeers; i++) {
                Ptr<PeerApp> peer = (*m_peerApps)[i];
                Simulator::Schedule(Seconds(0.01 + (0.001 * (double)i)), &PeerApp::GossipSignature, peer, peer->GetPeerId());
            }
        }
    }

    void CheckForNextRound() {
        if (m_globalAgg) {
            uint32_t globalCurrentRound = m_globalAgg->GetCurrentRound();

            // If global aggregator moved to a new round, update our state
            if (globalCurrentRound > m_currentRound) {
                StartNextRound(globalCurrentRound);

                // Broadcast round change to worker processes
                if (m_mpiMaster) {
                    BroadcastRoundChange(globalCurrentRound);
                }
            }
        }

        // Schedule next check
        Simulator::Schedule(Seconds(0.1), &SubnetAggregatorApp::CheckForNextRound, this);
    }

    void BroadcastRoundChange(uint32_t newRound) {
        // Only the master process should broadcast
        if (!m_mpiMaster)
            return;

        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        // Broadcast to all worker processes
        for (int i = 1; i < world_size; i++) {
            MPI_Send(&newRound, 1, MPI_UINT32_T, i, 999, MPI_COMM_WORLD);
        }

        LogMPI("Master process broadcast round change to " + std::to_string(newRound), 0);
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
    uint32_t m_currentRound;
    uint32_t m_totalRounds;
    std::set<uint32_t> m_receivedPeers;
    Ptr<GlobalAggregatorApp> m_globalAgg;
    std::vector<Ptr<PeerApp> > *m_peerApps; // Pointer to all peer apps in the subnet
    int m_mpiRank;
    bool m_mpiMaster;
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

            // Use the new logging system with process information
            LogMPI("Master received subnet aggregation via MPI from subnet " +
                  std::to_string(data), status.MPI_SOURCE);
        }
    }

    // Schedule next check
    Simulator::Schedule(Seconds(0.01), &CheckForMPIMessages, globalAggApp, world_size);
}

// Receive messages from master process (in worker processes)
void CheckForMPIRoundChanges(std::vector<Ptr<SubnetAggregatorApp>> &subnetAggApps, int world_rank) {
    // Only worker processes should check for round change messages
    if (world_rank == 0)
        return;

    MPI_Status status;
    int flag = 0;
    uint32_t newRound;

    // Check for round change messages from the master process
    MPI_Iprobe(0, 999, MPI_COMM_WORLD, &flag, &status);

    if (flag) {
        // Receive the round change message
        MPI_Recv(&newRound, 1, MPI_UINT32_T, 0, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Use special process-specific logging
        LogMPI("Process " + std::to_string(world_rank) + " received round change to " +
              std::to_string(newRound), world_rank);

        // Notify all subnet aggregators on this process about the round change
        for (auto aggApp : subnetAggApps) {
            aggApp->StartNextRound(newRound);
        }
    }

    // Schedule next check
    Simulator::Schedule(Seconds(0.02), &CheckForMPIRoundChanges, std::ref(subnetAggApps), world_rank);
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
    uint32_t nSubnets = 5; // Default total number of subnets across all processes
    uint32_t nPeersPerSubnet = 512; // Default number of peers per subnet
    uint32_t nRounds = 5; // Default number of rounds of aggregation
    uint32_t fanOutConn = 8; // Number of P2P connections per peer (fan-out)

    // Parse command line arguments
    CommandLine cmd;
    cmd.AddValue("nSubnets", "Total number of subnets across all processes", nSubnets);
    cmd.AddValue("nPeersPerSubnet", "Number of peers per subnet", nPeersPerSubnet);
    cmd.AddValue("nRounds", "Number of rounds of aggregation to run sequentially", nRounds);
    cmd.AddValue("fanOutConn", "Number of P2P connections per peer", fanOutConn);
    cmd.Parse(argc, argv);

    // Calculate which subnets to process on this rank
    uint32_t subnets_per_rank = nSubnets / world_size;
    uint32_t extra_subnets = nSubnets % world_size;
    uint32_t my_subnet_start = world_rank * subnets_per_rank + std::min(world_rank, static_cast<int>(extra_subnets));
    uint32_t my_subnet_count = subnets_per_rank + (world_rank < static_cast<int>(extra_subnets) ? 1 : 0);
    uint32_t my_subnet_end = my_subnet_start + my_subnet_count;

    if (world_rank == 0) {
        std::cout << "Running simulation with " << world_size << " MPI processes" << std::endl;
        std::cout << "Total subnets: " << nSubnets << ", peers per subnet: " << nPeersPerSubnet << ", rounds: " << nRounds << std::endl;
        std::cout << "Each peer will have up to " << fanOutConn << " P2P connections" << std::endl;
    }

    std::cout << "Process " << world_rank << " handling subnets " << my_subnet_start
              << " to " << (my_subnet_end-1) << " (" << my_subnet_count << " subnets)" << std::endl;

    // Global aggregator is needed only on the master process (rank 0)
    Ptr<GlobalAggregatorApp> globalAggApp;
    NodeContainer globalAggNode;

    if (world_rank == 0) {
        globalAggNode.Create(1);
        globalAggApp = CreateObject<GlobalAggregatorApp>();
        globalAggApp->Setup(nSubnets, nRounds);
        globalAggNode.Get(0)->AddApplication(globalAggApp);

        // Schedule regular checks for MPI messages
        Simulator::Schedule(Seconds(0.1), &CheckForMPIMessages, globalAggApp, world_size);
    }

    // Each process creates only the subnets it's responsible for
    std::vector<Ptr<SubnetAggregatorApp>> subnetAggApps(my_subnet_count);
    std::vector<NodeContainer> subnetNodes(my_subnet_count);
    std::vector<std::vector<Ptr<PeerApp>>> peerApps(my_subnet_count);
    std::vector<InternetStackHelper> internetStacks(my_subnet_count);
    std::vector<std::vector<Ipv4Address>> peerAddresses(my_subnet_count);

    // Create only subnets assigned to this process
    for (uint32_t i = 0; i < my_subnet_count; i++) {
        uint32_t s = my_subnet_start + i;  // Global subnet ID
        subnetNodes[i].Create(nPeersPerSubnet + 1);  // +1 for aggregator

        // Install Internet stack on all nodes in this subnet
        internetStacks[i].Install(subnetNodes[i]);

        // Set up a subnet-specific addressing scheme to avoid IP collisions
        Ipv4AddressHelper address;
        std::ostringstream subnetAddr;
        // Create unique subnet addresss (10.s.0.0)
        subnetAddr << "10." << s << ".0.0";
        address.SetBase(subnetAddr.str().c_str(), "255.255.0.0");

        // Use a point to point helper for this subnet
        PointToPointHelper p2p;
        p2p.SetDeviceAttribute("DataRate", StringValue("5Mbps"));
        p2p.SetChannelAttribute("Delay", StringValue("2ms"));

        // Create storage for peer addresses and SignatureServer instances
        peerAddresses[i].resize(nPeersPerSubnet);
        std::vector<Ptr<SignatureServer>> signatureServers(nPeersPerSubnet);

        // First, set up all nodes with internet and assign addresses
        for (uint32_t p = 0; p < nPeersPerSubnet; p++) {
            // Create a "stub" interface for each peer
            NetDeviceContainer stubDevices = p2p.Install(NodeContainer(subnetNodes[i].Get(p), subnetNodes[i].Get(nPeersPerSubnet)));
            Ipv4InterfaceContainer interfaces = address.Assign(stubDevices);
            peerAddresses[i][p] = interfaces.GetAddress(0);

            // Create a signature server for each peer
            signatureServers[p] = CreateObject<SignatureServer>();
            subnetNodes[i].Get(p)->AddApplication(signatureServers[p]);
        }

        // Create the subnet aggregator app
        Ptr<SubnetAggregatorApp> aggApp = CreateObject<SubnetAggregatorApp>();
        aggApp->Setup(s, nPeersPerSubnet, globalAggApp, &peerApps[i], world_rank);
        subnetNodes[i].Get(nPeersPerSubnet)->AddApplication(aggApp);
        subnetAggApps[i] = aggApp;

        // Create peer apps for all nodes and set up P2P connections
        peerApps[i].resize(nPeersPerSubnet);
        for (uint32_t p = 0; p < nPeersPerSubnet; p++) {
            Ptr<PeerApp> peerApp = CreateObject<PeerApp>();
            peerApp->Setup(p, aggApp, &peerApps[i], nPeersPerSubnet);
            peerApp->SetSignatureServer(signatureServers[p]);
            subnetNodes[i].Get(p)->AddApplication(peerApp);
            peerApps[i][p] = peerApp;
        }

        // Now set up the P2P connections between peers using a "smart" connection strategy:
        // - Each peer connects to fanOutConn other peers
        // - We use a structured topology to ensure good network coverage
        // - Avoid creating O(n²) connections but maintain good message propagation
        for (uint32_t p = 0; p < nPeersPerSubnet; p++) {
            std::set<uint32_t> connections;

            // 1. Connect to the next fanOutConn/2 sequential peers (with wraparound)
            for (uint32_t j = 1; j <= fanOutConn/2; j++) {
                uint32_t targetId = (p + j) % nPeersPerSubnet;
                connections.insert(targetId);
            }

            // 2. Connect to fanOutConn/2 exponentially distributed peers for "small world" properties
            for (uint32_t j = 0; j < fanOutConn/2; j++) {
                // Use a simple power-law distribution to select long-range peers
                uint32_t jump = (1 << (j + 1)) % nPeersPerSubnet; // 2, 4, 8, 16, ...
                uint32_t targetId = (p + jump) % nPeersPerSubnet;
                connections.insert(targetId);

                // If we got a duplicate (already in connections), try one more random peer
                if (connections.size() <= fanOutConn/2 + j) {
                    uint32_t randomPeer = rand() % nPeersPerSubnet;
                    if (randomPeer != p) {
                        connections.insert(randomPeer);
                    }
                }
            }

            // Now create the actual connections for this peer
            for (uint32_t targetId : connections) {
                if (targetId != p) { // Don't connect to self
                    peerApps[i][p]->AddPeerConnection(targetId, peerAddresses[i][targetId]);
                }
            }
        }
    }

    // Schedule regular checks for MPI round changes in worker processes
    if (world_rank != 0) {
        Simulator::Schedule(Seconds(0.1), &CheckForMPIRoundChanges, std::ref(subnetAggApps), world_rank);
    }

    // Set longer timeout to ensure simulation completes all rounds
    Simulator::Stop(Seconds(60.0));  // Increased from 10.0 to 60.0 seconds

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

    // Use MPI_Reduce to collect total messages sent across all processes
    uint64_t local_messages_sent = PeerApp::GetTotalMessagesSent();
    uint64_t total_messages_sent = 0;
    MPI_Reduce(&local_messages_sent, &total_messages_sent, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    // Print statistics on the master process
    if (world_rank == 0) {
        std::cout << "\n--- SIMULATION STATISTICS ---" << std::endl;
        std::cout << "MPI processes: " << world_size << std::endl;
        std::cout << "Total virtual time: " << std::fixed << std::setprecision(3) << max_virtualTimeSeconds << " seconds" << std::endl;
        std::cout << "Total real execution time: " << std::fixed << std::setprecision(3) << max_realTimeSeconds << " seconds" << std::endl;
        std::cout << "Ratio (virtual/real): " << std::fixed << std::setprecision(6) << (max_virtualTimeSeconds / max_realTimeSeconds) << std::endl;
        std::cout << "Total number of peers: " << nSubnets * nPeersPerSubnet << std::endl;
        std::cout << "Average connections per peer: " << fanOutConn << std::endl;
        std::cout << "Total messages sent: " << total_messages_sent << std::endl;
        std::cout << "----------------------------" << std::endl;
    }

    Simulator::Destroy();

    // Clean up MPI datatype
    MPI_Type_free(&mpi_signature_type);

    MPI_Finalize();  // Finalize MPI
    return 0;
}

// Signature Server implementation
TypeId
SignatureServer::GetTypeId(void)
{
  static TypeId tid = TypeId("ns3::SignatureServer")
    .SetParent<Application>()
    .AddConstructor<SignatureServer>();
  return tid;
}

SignatureServer::SignatureServer()
  : m_socket(0),
    m_port(9000)
{
}

SignatureServer::~SignatureServer()
{
}

void
SignatureServer::DoDispose(void)
{
  m_socket = 0;
  Application::DoDispose();
}

void
SignatureServer::StartApplication(void)
{
  if (m_socket == nullptr) {
    TypeId tid = TypeId::LookupByName("ns3::UdpSocketFactory");
    m_socket = Socket::CreateSocket(GetNode(), tid);
    InetSocketAddress local = InetSocketAddress(Ipv4Address::GetAny(), m_port);
    m_socket->Bind(local);
    m_socket->SetRecvCallback(MakeCallback(&SignatureServer::HandleRead, this));
  }
}

void
SignatureServer::StopApplication(void)
{
  if (m_socket != nullptr) {
    m_socket->Close();
    m_socket->SetRecvCallback(MakeNullCallback<void, Ptr<Socket>>());
    m_socket = 0;
  }
}

void
SignatureServer::SetPeerApp(Ptr<PeerApp> peerApp)
{
  m_peerApp = peerApp;
}

void
SignatureServer::HandleRead(Ptr<Socket> socket)
{
  Ptr<Packet> packet;
  Address from;

  while (packet = socket->RecvFrom(from)) {
    // Extract signature header
    SignatureHeader sigHeader;
    packet->RemoveHeader(sigHeader);

    // Create a SignatureMessage to be compatible with existing code
    SignatureMessage sig(sigHeader.GetOriginPeerId(), 0);

    // Process the received signature via the parent PeerApp
    if (m_peerApp) {
      m_peerApp->ReceiveGossipSignature(sig);
    }
  }
}

// Signature Client implementation
TypeId
SignatureClient::GetTypeId(void)
{
  static TypeId tid = TypeId("ns3::SignatureClient")
    .SetParent<Application>()
    .AddConstructor<SignatureClient>();
  return tid;
}

SignatureClient::SignatureClient()
  : m_socket(0),
    m_peerPort(9000),
    m_connected(false),
    m_lastOriginPeerId(0)
{
}

SignatureClient::~SignatureClient()
{
}

void
SignatureClient::DoDispose(void)
{
  m_socket = 0;
  Application::DoDispose();
}

void
SignatureClient::SetAttribute(std::string name, const AttributeValue &value)
{
  if (name == "PeerAddress") {
    const AddressValue *addressValue = dynamic_cast<const AddressValue*>(&value);
    if (addressValue) {
      m_peerAddress = addressValue->Get();
    }
  } else if (name == "PeerPort") {
    const UintegerValue *portValue = dynamic_cast<const UintegerValue*>(&value);
    if (portValue) {
      m_peerPort = portValue->Get();
    }
  }
}

void
SignatureClient::Setup(Ipv4Address address, uint16_t port)
{
  m_peerAddress = InetSocketAddress(address, port);
  m_peerPort = port;
}

void
SignatureClient::StartApplication(void)
{
  // Create socket and connect to peer
  if (m_socket == nullptr) {
    TypeId tid = TypeId::LookupByName("ns3::UdpSocketFactory");
    m_socket = Socket::CreateSocket(GetNode(), tid);

    if (Ipv4Address::IsMatchingType(m_peerAddress)) {
      m_socket->Connect(InetSocketAddress(Ipv4Address::ConvertFrom(m_peerAddress), m_peerPort));
    }

    m_socket->SetConnectCallback(
      MakeCallback(&SignatureClient::ConnectionSucceeded, this),
      MakeCallback(&SignatureClient::ConnectionFailed, this));

    m_connected = true;
  }
}

void
SignatureClient::StopApplication(void)
{
  if (m_socket == nullptr) {
    m_socket->Close();
    m_socket = 0;
    m_connected = false;
  }
}

void
SignatureClient::ConnectionSucceeded(Ptr<Socket> socket)
{
  m_connected = true;
  // If there was a pending message, try to send it now
  if (m_lastOriginPeerId > 0) {
    SendSignature(m_lastOriginPeerId);
    m_lastOriginPeerId = 0;
  }
}

void
SignatureClient::ConnectionFailed(Ptr<Socket> socket)
{
  m_connected = false;
}

void
SignatureClient::SendSignature(uint32_t originPeerId)
{
  if (!m_connected) {
    m_lastOriginPeerId = originPeerId;
    return;
  }

  Ptr<Packet> packet = Create<Packet>();
  SignatureHeader header;
  header.SetOriginPeerId(originPeerId);

  packet->AddHeader(header);

  PeerApp::s_totalMessagesSent++;
  // if (m_socket->Send(packet) >= 0) {
  // }
}

