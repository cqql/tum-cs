#!/bin/bash

set -e -u

OTHERN="${OTHERN:-zbox07}"
OCD="${OCD:-m}"
RATE="${RATE:-20}"
MCS="${MCS:-0}"
COUNT="${COUNT:-1000}"
LOG="${LOG:-log}"
FREQ="${FREQ:-5700M}"

cd
cd "$OCD"

./ptmbeacon -a de:ad:be:ef:42:01 -h HT20 -l0 -c"$COUNT" -b0 -r0 -i10.0.0.1 wlan0 "$FREQ" >"$LOG" &
bgp=$!

sleep .5

ssh "$OTHERN" "pkill lt-ptmbeacon; cd \"$OCD\"; sleep .3; ./ptmbeacon -a de:ad:be:ef:42:02 -h HT20 -l0 -c\"$COUNT\" -b\"$RATE\" -r\"$MCS\" -i10.0.0.2 wlan0 \"$FREQ\" >/dev/null"

sleep .1

kill -INT $bgp
wait $bgp

pc=$(grep -c "" "$LOG" | cut -d\  -f1)
echo "LQ: $pc / $COUNT - "$(bc <<<"scale=4; $pc / $COUNT")
