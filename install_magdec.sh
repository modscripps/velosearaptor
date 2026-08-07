#!/usr/bin/env bash
# Build the magdec executable used to calculate magnetic declination.
#
# Requires a C compiler and make. The resulting geomag/ directory is gitignored
# and is found by velosearaptor at runtime (see madcp._find_magdec).
set -euo pipefail

cd "$(dirname "$0")"

if [ -d geomag ]; then
    echo "geomag/ already exists - remove it first to rebuild."
    exit 1
fi

git clone https://currents.soest.hawaii.edu/git/Oceanography_Tools/geomag
cd geomag
# Drop the version control stuff, we only need the built executable.
rm -rf .git .gitignore
make magdec

echo "magdec built at $(pwd)/magdec"
