#!/bin/bash

d=$(date +%s)

for f in 5700 4920 5240; do
	for r in 2 100 20000; do
		for m in 2 3; do
			dur=$(( $r * 1200 ))
			FREQ=$f"M" COUNT=$dur RATE=$r LOG=maes-$d-$f-$r-$m MCS=$m ./domeasure.sh
		done
	done
done
