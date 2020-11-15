#!/bin/bash

cargo run --release
gnuplot -persist plot.gpt
