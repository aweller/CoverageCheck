import os

def unite(output, target_folder = "./"):
	files = [x for x in os.listdir(target_folder) if x.endswith("_coverage.tsv")]
	
	with open(target_folder + output, "w") as out:
		for filename in files:
			sample = filename[:-13]
		
			with open(target_folder + filename) as handle:
				for row in handle:
					result = sample +"\t"+ row 
					out.write(result)
	return target_folder + output
	
####################################################################################

if __name__ == '__main__':

	files = [x for x in os.listdir(".") if x.endswith("_coverage.tsv")]
	
	for filename in files:
	
		sample = filename[:-13]
	
		with open(filename) as handle:
			for row in handle:
				print sample +"\t"+ row, 
