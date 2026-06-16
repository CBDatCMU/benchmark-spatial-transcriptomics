#!/bin/bash
# Downloads example_spinal_chord_inactive.zip from Zenodo (record 11120922),
# extracts it into data/, and removes the zip.

# Skips download if data/example_spinal_chord_inactive/ already exists and is non-empty.

set -euo pipefail

ZENODO_URL="https://zenodo.org/records/11120922/files/example_spinal_chord_inactive.zip"
ZIP_NAME="example_spinal_chord_inactive.zip"
DATA_DIR="$(cd "$(dirname "$0")" && pwd)/data"
TARGET_DIR="${DATA_DIR}/example_spinal_chord_inactive"

echo "Target: ${TARGET_DIR}"

# Skip if already present
if [ -d "${TARGET_DIR}" ] && [ -n "$(ls -A "${TARGET_DIR}" 2>/dev/null)" ]; then
    echo "Data already present, skipping download."
    exit 0
fi

mkdir -p "${DATA_DIR}"

echo "Downloading ${ZIP_NAME} (1.7 GB)..."
curl -L --progress-bar -o "${DATA_DIR}/${ZIP_NAME}" "${ZENODO_URL}"

echo "Extracting..."
unzip -q "${DATA_DIR}/${ZIP_NAME}" -d "${DATA_DIR}"

rm "${DATA_DIR}/${ZIP_NAME}"
echo "Done. Data path: ${TARGET_DIR}"