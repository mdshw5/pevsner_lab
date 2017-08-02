#!/usr/bin/env bash

bfc -s 3g -t 2 <(seqtk mergepe $1 $2) <(seqtk mergepe $1 $2) > $3