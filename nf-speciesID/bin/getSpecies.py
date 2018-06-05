#!/usr/bin/env python3
import sys
import os
import shutil
# Expected file name: {pair_id}.centrifuge_report.tsv

input_file = sys.argv[1] #Input centrifuge report file
class_folder = sys.argv[2] #Input path to save taxonomized file

_prefix = input_file.split(".")[0]

with open(input_file, "r") as f:
    next(f)
    _name = []
    _reads = []
    for line in f:
        name, _, _, _, _, reads, _ = line.split("\t")
        _name.append(name)
        _reads.append(int(reads))

_max_index = _reads.index(max(_reads))
species_name = _name[_max_index]
with open("{}_speciation.csv".format(_prefix), "w+") as output:
    output.write("{},{}\n".format(_prefix, species_name))

# Move reads into folder
species_name = species_name.replace(" ", "_")


species_name_folder = os.path.join(class_folder, species_name)
if not os.path.isdir(species_name_folder):
    os.makedirs(species_name_folder)

gz_files = os.listdir(".")
for gz in gz_files:
    if gz.endswith(".gz"):
        shutil.move(gz, os.path.join(species_name_folder, os.path.basename(gz)))
