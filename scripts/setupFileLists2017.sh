#!/bin/env bash

pushd () {
    command pushd "$@" > /dev/null
}

popd () {
    command popd "$@" > /dev/null
}

printf "Updating filelists for datasets used ...\n"
printf "First deleting old filelists ...\n"

if [ -z "${TQZ_TOOLS_PATH}" ]; then
    printf "TQZ_TOOLS_PATH not set, setting to ${PWD}\n"
    TQZ_TOOLS_PATH="${PWD}"
fi

rm -rf "${TQZ_TOOLS_PATH}"/configs/2017/datasets/fileLists/*

printf "Done!\n"
printf "Now outputting the lists of the dataset directories into their relevant files ...\n"

# Nominal skims
pushd /scratch/data/tZqSkimsRun2017/
for i in *; do
    if [ -d "${i}" ]; then
        pushd "${i}"
        ls -d "${PWD}"/* > "${TQZ_TOOLS_PATH}/configs/2017/datasets/fileLists/${i}Files.txt"
        popd
    fi
done
popd

# Overwrite with post trigger skims where they exist
pushd /data0/data/TopPhysics/postTriggerSkims2017/
for i in *; do
    if [ -d "${i}" ]; then
        pushd "${i}"
        ls -d "${PWD}"/* > "${TQZ_TOOLS_PATH}/configs/2017/datasets/fileLists/${i}Files.txt"
        popd
    fi
done
popd

printf "Done!\n"
printf "Filelists have been updated.\n"
