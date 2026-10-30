#!/bin/sh
# Download the latest OpenNGC catalog (NGC, IC, Messier) into ./data
set -e
cd "$(dirname "$0")"
mkdir -p data
base=https://raw.githubusercontent.com/mattiaverga/OpenNGC/master/database_files
curl -sSfo data/NGC.csv "$base/NGC.csv"
curl -sSfo data/addendum.csv "$base/addendum.csv"
echo "Catalog updated: $(wc -l < data/NGC.csv) NGC/IC rows"
