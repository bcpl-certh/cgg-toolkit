#!/usr/bin/env bash

printf "Setting environmental PRSSDIR\n"
export PRSSDIR=$(pwd)
printf "PRSSDIR=%s\n" $(printenv PRSSDIR)

./generage "$@"
