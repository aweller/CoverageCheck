import os

files = [x for x in os.listdir(".") if x.endswith("_coverage.tsv")]

for filename in files:

	sample = filename[:-13]

	with open(filename) as handle:
		for row in handle:
			print sample +"\t"+ row, 
