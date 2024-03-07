#!/bin/bash

## set the directory of loco-pipe as PIPEDIR
cd `dirname ${0}`/..
PIPEDIR=`pwd -P`
echo $PIPEDIR

## change the timestamp of the index file of the reference genome
touch $PIPEDIR/toyfish/reference/toy_refgen.fa.fai 

## update the path in the config file
mkdir -p $PIPEDIR/toyfish/config
sed "s|/path/to/loco-pipe|$PIPEDIR|g" $PIPEDIR/config.yaml > $PIPEDIR/toyfish/config/config.yaml

## update the path in the sample table
sed -i "s|/path/to/loco-pipe|$PIPEDIR|g" $PIPEDIR/toyfish/docs/sample_table.tsv

