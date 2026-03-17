#!/bin/bash

rm -rf ./results/baseline/
python3 scripts/trace_exec_training_list.py --trace_dir ./traces --results_dir ./results/baseline
cp ./results/results.csv ReportGenerators/