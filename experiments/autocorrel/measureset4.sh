#!/bin/bash

d=$(date +%s)

for fo in "5240 11"; do
	for i in {1..100}; do
		m="${fo/* /}"
		f="${fo/ */}"
		FREQ=$f"M" COUNT=5000 RATE=5 BC=150 LOG=heaviside-$f-$m-$i MCS=$m ./domeasure.sh
	done
done
