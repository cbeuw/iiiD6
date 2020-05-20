#!/bin/bash
cargo build --release;
mkdir -p render;
cd render;
for n in {2..4}
do
    for l in $(seq 0 $(($n-1)))
    do
        for m in $(seq 0 $l)
        do
            echo "Rendering n=$n, l=$l, m=$m"
            ../target/release/iiiD6 $n $l $m
        done
    done
done