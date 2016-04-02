#!/usr/bin/python

import sys

PATTERN='%%REPORT_ID%%'
REPORT_ID=""

if len(sys.argv)>1:
  REPORT_ID=sys.argv[1]

for line in sys.stdin:
  line = line.rstrip("\n")
  if PATTERN in line:
    line = line.replace(PATTERN, REPORT_ID)
  print line

