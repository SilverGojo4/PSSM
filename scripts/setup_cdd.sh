#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# ===================== Usage =====================
if [ $# -lt 1 ]; then
  echo "Usage: $0 <PROJECT_DIR>"
  exit 1
fi

PROJECT_DIR=$1
CDD_TAR="$PROJECT_DIR/blastdb/cdd.tar"
CDD_DIR="$PROJECT_DIR/blastdb/cdd"

# ===================== Step 1: Check Files =====================
if [ ! -f "$CDD_TAR" ]; then
  echo "❌ Missing $CDD_TAR"
  echo "Please download it first from:"
  echo "  ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cdd.tar.gz"
  exit 1
fi

mkdir -p "$CDD_DIR"

# ===================== Step 2: Extract Database =====================
echo "📦 Extracting CDD from $CDD_TAR ..."
tar -xvf "$CDD_TAR" -C "$CDD_DIR"
echo "✅ Extraction complete."

# ===================== Step 3: Build RPS-BLAST Database =====================
echo "🔨 Building RPS-BLAST database ..."
cd "$CDD_DIR"

# 確認是否存在 .smp 檔
if ! ls *.smp 1> /dev/null 2>&1; then
  echo "❌ No .smp files found in $CDD_DIR"
  exit 1
fi

# 自動建立 profile list (Cdd.pn)
ls *.smp > Cdd.pn

# 使用 makeprofiledb 建立 RPS-BLAST 資料庫
makeprofiledb -in Cdd.pn -out Cdd -dbtype rps -title "Conserved Domain Database (CDD)"

# ===================== Step 4: Verify =====================
echo "📂 Verifying generated files ..."
ls -lh Cdd.*

echo "✅ CDD setup completed successfully!"
