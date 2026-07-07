#!/bin/bash
# Downloads four mouse brain Xenium datasets from Zenodo (record 10566172).
# Skips download if target directory already exists and is non-empty.
set -euo pipefail

BASE_URL="https://zenodo.org/records/10566172/files"
DATA_DIR="$(cd "$(dirname "$0")" && pwd)/data"

download_dataset() {
    local ZIP_NAME="$1"
    local SIZE="$2"
    local TARGET_DIR="${DATA_DIR}/${ZIP_NAME%.zip}"

    if [ -d "${TARGET_DIR}" ] && [ -n "$(ls -A "${TARGET_DIR}" 2>/dev/null)" ]; then
        echo "${ZIP_NAME%.zip}: already present, skipping."
        return
    fi

    mkdir -p "${DATA_DIR}"
    echo "Downloading ${ZIP_NAME} (${SIZE})..."
    curl -L --progress-bar -o "${DATA_DIR}/${ZIP_NAME}" "${BASE_URL}/${ZIP_NAME}"
    echo "Extracting..."
    unzip -q "${DATA_DIR}/${ZIP_NAME}" -d "${DATA_DIR}"
    rm "${DATA_DIR}/${ZIP_NAME}"
    echo "Done. Data path: ${TARGET_DIR}"
}

download_dataset "output-XETG00047__0003467__brain1__20230331__153000.zip" "6.1 GB"
download_dataset "output-XETG00047__0003467__brain2__20230331__153000.zip" "6.8 GB"
download_dataset "output-XETG00047__0003304__brain3__20230331__153000.zip" "6.5 GB"
download_dataset "output-XETG00047__0003304__brain4__20230331__153000.zip" "7.0 GB"
