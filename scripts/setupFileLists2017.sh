#!/bin/env bash

set -Eeuo pipefail
shopt -s extglob

pushd () {
    command pushd "$@" > /dev/null
}

popd () {
    command popd "$@" > /dev/null
}

TQZ_TOOLS_PATH="${TZQ_TOOLS_PATH-$PWD}"

printf "Updating filelists for datasets used ...\n"
printf "First deleting old filelists ...\n"

rm -rf "${TQZ_TOOLS_PATH}"/configs/2017/datasets/fileLists/*

printf "Done!\n"
printf "Now outputting the lists of the dataset directories into their relevant files ...\n"

# Nominal skims
pushd /scratch/data/tZqSkimsRun2017/
for i in *; do
    if [ -d "${i}" ]; then
        pushd "${i}"
        ls -d "${PWD}"/* >> "${TQZ_TOOLS_PATH}/configs/2017/datasets/fileLists/${i%_ext+([[:digit:]])}Files.txt" || printf "WARNING: ${PWD} is empty\n"
        popd
    fi
done
popd

# Overwrite with post trigger skims where they exist
pushd /data0/data/TopPhysics/postTriggerSkims2017/
for i in *; do
    if [ -d "${i}" ]; then
        pushd "${i}"
        ls -d "${PWD}"/* > "${TQZ_TOOLS_PATH}/configs/2017/datasets/fileLists/${i}Files.txt" || printf "WARNING: ${PWD} is empty\n"
        popd
    fi
done
popd

printf "Done!\n"
printf "Filelists have been updated.\n"
