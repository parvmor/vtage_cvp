#!/usr/bin/env python3

import os
import subprocess
import sys

ipc = 1.0
trace_cnt = 0

if len(sys.argv) < 2:
    print('usage: ./get_ipc.py <dir> [<traces-dir>]')
    sys.exit(1)
elif len(sys.argv) == 2:
    sys.argv.append('./traces')
traces = os.listdir(sys.argv[2])

for tr in traces:
    print('running: ' + tr + '...')
    cmd = './' + sys.argv[1] + '/cvp'
    opt = '-F 16,0,0,0,0 -f 5 -M 0 -A 0 -w 256 -v ' + sys.argv[2] + '/' + tr
    output = subprocess.Popen(cmd + ' ' + opt, shell=True, stdout=subprocess.PIPE).stdout.read()
    cur_ipc = float(output.decode('utf-8').split('\n')[-7].split('=')[1].strip())
    ipc, trace_cnt = ipc * cur_ipc, trace_cnt + 1

print(ipc ** (1.0 / trace_cnt))
