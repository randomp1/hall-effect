#!/bin/bash

# Save data, generate all graphs
./saveResults.py
./analysis.py
./analysis2.py

# Delete unnecessary graphs from graph folders
cd ./graphs_raw
./rm_unnecessary.sh
cd ../correctedGraphs
./rm_unnecessary.sh
cd ../fitnessGraphs
./rm_unnecessary.sh

cd ..

# Remove old archive if present, zip all graphs into single archive
rm ./allGraphs.zip
zip allGraphs ./graphs_raw/*.pdf ./correctedGraphs/*.pdf ./fitnessGraphs/*.pdf
