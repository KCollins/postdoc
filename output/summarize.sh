#!/bin/sh
for y in $(seq -w 2012 2018); do
for m in $(seq -w 01 12); do
for d in $(seq -w 01 31); do
convert +append ${y}-${m}-${d}*H.png ${y}-${m}-${d}*E.png ${y}-${m}-${d}-summary.png;
done;
done;
done

for y in $(seq -w 2012 2018); do
for m in $(seq -w 01 12); do
for d in $(seq -w 01 31); do
convert +append ${y}-${m}-${d}-summary.png ${y}-${m}-${d}*Z.png ${y}-${m}-${d}-summary.png;
done; done; done
