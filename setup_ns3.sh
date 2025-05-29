#!/bin/bash
set -e  # Exit on error

echo "Setting up NS-3 for the Gossip Protocol Simulation"
echo "===================================================="

# Create directory for external dependencies
mkdir -p external
cd external

# Set NS-3 version
NS3_VERSION="3.44"
NS3_DIR="ns-allinone-$NS3_VERSION"
NS3_ARCHIVE="$NS3_DIR.tar.bz2"
NS3_URL="https://www.nsnam.org/releases/$NS3_ARCHIVE"

# Download NS-3 if not already downloaded
if [ ! -f "$NS3_ARCHIVE" ]; then
    echo "Downloading NS-3 version $NS3_VERSION..."
    curl -O "$NS3_URL" || wget "$NS3_URL"
else
    echo "NS-3 archive already downloaded."
fi

# Extract NS-3 if not already extracted
if [ ! -d "$NS3_DIR" ]; then
    echo "Extracting NS-3..."
    tar -xf "$NS3_ARCHIVE"
else
    echo "NS-3 already extracted."
fi

# Enter NS-3 directory and build
cd "$NS3_DIR"
if [ ! -d "install" ]; then
    echo "Building NS-3 (this may take a while)..."

    # Configure NS-3 with CMake
    mkdir -p build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DENABLE_EXAMPLES=OFF \
          -DENABLE_TESTS=OFF \
          -DNS3_BINDINGS_INSTALL=OFF \
          -DCMAKE_INSTALL_PREFIX=../install \
          ../ns-$NS3_VERSION

    # Build and install NS-3
    make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2) install

    echo "NS-3 built successfully!"
else
    echo "NS-3 already built."
fi

cd "../../.."  # Return to project root directory

echo ""
echo "NS-3 setup complete!"
echo "You can now build the project with: mkdir -p build && cd build && cmake .. && make"
