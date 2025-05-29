# Gossip Protocol Simulation in NS-3

## Overview

This project implements a distributed gossip protocol for signature distribution in a hierarchical network using the NS-3 network simulation framework. The simulation is parallelized using MPI (Message Passing Interface) to distribute computations across multiple processes, allowing for efficient simulation of large-scale networks.

## System Architecture

The system is organized as a hierarchical network with three levels:

1. **Peers**: Basic entities that create and gossip signatures.
2. **Subnet Aggregators**: Collect signatures from peers within their subnet and forward them to the global aggregator.
3. **Global Aggregator**: Receives aggregations from all subnet aggregators.

### Gossip Protocol

Instead of peers sending signatures directly to subnet aggregators, signatures are gossiped between peers with a fan-out factor of 6. Each peer:
- Generates its own signature
- Forwards new signatures to 6 randomly selected peers in the same subnet
- Only forwards each unique signature once

Subnet aggregators monitor the signatures spreading through the network and collect them until reaching a threshold (2/3 + 1 of all peers in the subnet).

## MPI Parallelization

The simulation leverages MPI to distribute the computational workload:

- Each MPI process handles a subset of subnets
- Only rank 0 (master process) hosts the global aggregator
- Workers communicate with the master through MPI messages
- Each process simulates only its assigned subnets, reducing memory usage and improving performance

## Requirements

- C++23 compiler
- CMake 3.25 or higher
- MPI implementation (OpenMPI recommended)

## Setup and Installation

### Step 1: Download and Build NS-3

The project includes a script to download, build, and install NS-3.44 locally within the project:

```bash
# Make the setup script executable (if not already)
chmod +x setup_ns3.sh

# Run the setup script
./setup_ns3.sh
```

This script will:
1. Create an `external` directory
2. Download NS-3.44 source code
3. Extract and build NS-3 with optimized settings
4. Install NS-3 within the project folder

The process may take 20-30 minutes depending on your system.

### Step 2: Build the Project

```bash
mkdir -p cmake-build-debug
cd cmake-build-debug
cmake ..
make
```

## Running

To run the simulation with MPI:

```bash
mpirun -np <num_processes> ./main [options]
```

Replace `<num_processes>` with the number of MPI processes you want to use.

### Command Line Parameters

You can configure the simulation using the following command-line options:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--nSubnets=<value>` | Number of subnets in the simulation | 5       |
| `--nPeersPerSubnet=<value>` | Number of peers in each subnet | 512     |

### Examples

```bash
# Run with 4 MPI processes using default settings
mpirun -np 4 ./main

# Run with 8 MPI processes, 50 subnets, 512 peers per subnet
mpirun -np 8 ./main --nSubnets=5 --nPeersPerSubnet=512

# Run with 16 MPI processes, 200 subnets, 2048 peers per subnet
mpirun -np 16 ./main --nSubnets=16 --nPeersPerSubnet=1024
```

## Configuration Parameters

The main configuration parameters are in `main.cpp`:

- `nSubnets`: Total number of subnets in the simulation
- `nPeersPerSubnet`: Number of peers in each subnet
- `m_fanOut`: Number of peers that each peer gossips to (default: 6)
- `m_threshold`: Threshold for subnet aggregators (default: 2/3 * nPeers + 1)

## Performance Metrics

The simulation outputs several performance metrics:

- Total virtual time: Time taken in the simulated environment
- Total real execution time: Actual wall-clock time
- Ratio (virtual/real): Efficiency of simulation
- Total number of peers: Size of the simulated network

## Logs Format

All logs include timestamps in milliseconds and identify the source:
- Peer messages show the peer ID
- Subnet Aggregator (SA) messages show the subnet ID
- Global Aggregator messages show when subnet aggregations are received
