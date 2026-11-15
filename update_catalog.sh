#!/bin/sh
# Download the latest catalogs into ./data:
#   OpenNGC (NGC, IC, Messier) and the HYG star database (constellation-line stars only)
set -e
cd "$(dirname "$0")"
mkdir -p data
base=https://raw.githubusercontent.com/mattiaverga/OpenNGC/master/database_files
curl -sSfo data/NGC.csv "$base/NGC.csv"
curl -sSfo data/addendum.csv "$base/addendum.csv"
echo "Catalog updated: $(wc -l < data/NGC.csv) NGC/IC rows"

# Stars used by the constellation lines, from HYG v4.1 (CC BY-SA 4.0).
# data/constellation_lines.fab and constellation_names.fab come from Stellarium (skycultures/western).
tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT
curl -sSfLo "$tmp" https://raw.githubusercontent.com/astronexus/HYG-Database/main/hyg/CURRENT/hygdata_v41.csv
python3 - "$tmp" <<'EOF'
import csv, sys
line_hips = set()
with open("data/constellation_lines.fab") as f:
    for ln in f:
        p = ln.split()
        if len(p) > 2:
            line_hips.update(int(h) for h in p[2:])
out = []
with open(sys.argv[1], newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if not r["hip"] or r["proper"] == "Sol":
            continue
        hip, mag = int(r["hip"]), float(r["mag"])
        if hip in line_hips:
            out.append([hip, f"{float(r['ra']) * 15:.5f}", f"{float(r['dec']):.5f}", f"{mag:.2f}",
                        r["proper"], r["bayer"], r["flam"], r["con"]])
out.sort(key=lambda x: float(x[3]))
with open("data/constellation_stars.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["hip", "ra", "dec", "mag", "name", "bayer", "flam", "con"])
    w.writerows(out)
print(f"Constellation stars updated: {len(out)} stars")
EOF
