#!/usr/bin/env bash

bfc -1s 3g -k 61 -t 2 $1 | gzip -1 > $2