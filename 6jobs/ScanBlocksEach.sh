#!/bin/bash
# Block scan with a fixed thickness per block: ScanBlocks.sh with EACH=1 (total CE thickness = count x thickness, folder <folder>_<thickness>mEach/<E>GeV/<n>block).
# Usage: bash ScanBlocksEach.sh <template macro in build/> <thickness of each block in m> [dz in m, default 0.025], e.g. bash ScanBlocksEach.sh run_0.mac 0.4
EACH=1 exec bash "$(dirname "$0")/ScanBlocks.sh" "$@"
