#!/usr/bin/env bash
set -ue

CONFIG_FILE=$1
OUT_DIR=$2

mkdir -p ${OUT_DIR}

# Run NOVOPlasty assembly
NOVOPlasty.pl -c ${CONFIG_FILE} -o ${OUT_DIR}

echo "NOVOPlasty assembly finished. Plastome sequences are in ${OUT_DIR}"
