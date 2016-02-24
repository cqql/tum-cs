#!/bin/bash

d=$(date +%s)

m=0
for f in 5700 4920 5240; do
	for r in "1e9 10000000" "100 240000" "2 3600"; do
		for l in 1000 500 0; do
			ct="${r/* /}"
			rr="${r/ */}"
			FREQ=$f"M" COUNT=$ct RATE=$rr LOG=maes-$d-$f-$rr-$m-$l MCS=$m ./domeasure.sh
		done
	done
done
