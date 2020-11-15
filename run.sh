#!/bin/bash

rm render.csv
cargo run --release
gnuplot -persist plot.gpt
