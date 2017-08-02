#!/usr/bin/env bash

bcftools convert -O v $1 | bgzip > $2