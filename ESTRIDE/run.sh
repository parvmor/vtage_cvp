#!/usr/bin/env bash

set -e

make
./cvp -F 16,0,0,0,0 -f 5 -M 0 -A 0 -w 256 -v ../traces/srv_42.gz
