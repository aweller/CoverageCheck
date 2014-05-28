# contains the Classes for CoverageCheck

import numpy as np

#####################################################################################################################

class BadRegion():
    
    """
    Each instance corresponds to one contiguos run of bases 
    """
    
    def __init__(self, row):
        f = row.strip().split("\t")
        self.chrom, self.start, self.stop, self.amplicon, na, self.strand,  self.pos, self.value = f
        self.gene = self.amplicon.split("_")[0]
        
        self.first_start = int(self.start)
        self.first_pos = int(self.pos)
        
        self.last_start = int(self.start)
        self.last_stop = int(self.stop)
        self.last_pos = int(self.pos)
        
        self.rows = []
        self.values = []
        self.rows.append(row)
        self.values.append(self.value)
         
    def add_row(self,row):
        f = row.strip().split("\t")
        self.chrom, self.start, self.stop, self.amplicon, na, self.strand,  self.pos, self.value = f
        
        self.last_start = int(self.start)
        self.last_stop = int( self.stop)
        self.last_pos = int( self.pos)
        
        self.rows.append(row)
        self.values.append(self.value)
        
    def print_output(self):
        
        self.region_start = self.first_start + self.first_pos
        self.region_stop = self.last_start + self.last_pos
        
        mean_value = round(np.mean([float(x) for x in self.values]), 2)
        region_length = len(self.values)
        
        output = [self.gene, self.chrom, self.region_start, self.region_stop, mean_value, region_length]
        
        return "\t".join([str(x) for x in output])

#####################################################################################################################

class ExpectedVariants():
    """
    Each instance corresponds to a list of expected variants
    """
    
    def __init__(self, filename, folder = "./", dp_cutoff=None, strandbias=None, samplename=None):
        self.dict = {}
        self.chromposes = []
        self.filename = filename
        self.sample = samplename
        self.folder = folder
        self.min_dp = dp_cutoff
        self.max_strand_ratio = strandbias
        
        self.gene_locs = {}
        self.gene_variants = {}
        self.gene_variants["no_gene"] = {"total":0, "discoverable":0}

        self.all_genes = []
        
        ###########################################################
        # parse the gene locations
        
        with open("/home/andreas/bioinfo/projects/nhs_coverage_tool/scripts/human_gene_locations.txt") as handle:
            for row in handle:
                f = row.strip().split("\t")
                chrom, start, stop, gene = f[0], int(f[1]), int(f[2]), f[3]
                if not self.gene_locs.get(chrom):
                    self.gene_locs[chrom] = []
                
                self.gene_locs[chrom].append([start, stop, gene])
                self.gene_variants[gene] = {"total":0, "discoverable":0}
                self.all_genes.append(gene)
        
        ###########################################################
        # add the variants
        
        self.header_found = False
        self.header_row = []
        
        with open(filename) as handle:
            for row in handle:
                if row[0] == "#": continue
                f = row.strip().split("\t")
                chrom, pos = f[:2]

                try:
                    pos = int(pos)
                except:
                    if not self.header_found:
                        self.header_row = f
                        self.header_found = True
                        continue

                    else:
                        print "Skipping variant row:  " + row,
                        continue
                
                chrom, pos = f[:2]
                chrompos = chrom + "\t" + pos        
                
                self.dict[chrompos] = {}
                self.dict[chrompos]["row"] = str(row)[:-1]        
                self.dict[chrompos]["dp"] = 0
                self.dict[chrompos]["sb"] = None        
                self.chromposes.append(chrompos)
        
                ##################################################
                
                self.dict[chrompos]["gene"] = "no_gene"
                for [start, stop, gene] in self.gene_locs[chrom]:
                    if start < int(pos) < stop:
                        self.dict[chrompos]["gene"] = gene
                        break            
        
    def add_coverage(self, chrompos, dp):
        self.dict[chrompos]["dp"] = dp

    def add_strand_ratio(self, chrompos, strand_ratio):
        self.dict[chrompos]["sb"] = strand_ratio        
    
    def print_output(self):
        
        outname = self.folder + self.sample + "_expected_variant_coverage.tsv"
        
        with open(outname, "w") as out:
            
            self.header_row.extend(["gene", "coverage", "coverage_ok", "strand_ratio", "strand_ratio_ok"])
            header = "\t".join([str(x) for x in self.header_row])
            out.write(header + "\n")
            
            for chrompos in self.chromposes:
                var = self.dict.get(chrompos)
                
                gene = self.dict[chrompos]["gene"]
                dp_pass = var["dp"] > self.min_dp
                
                strand_ratio_pass = False
                if var.get("sb"):
                    strand_ratio_pass = var["sb"] < self.max_strand_ratio
                
                result = [var["row"], gene, var["dp"], dp_pass, var["sb"], strand_ratio_pass]
                result = "\t".join([str(x) for x in result])
                
                out.write(result + "\n")

                self.gene_variants[gene]["total"] += 1
                if dp_pass and strand_ratio_pass:
                    self.gene_variants[gene]["discoverable"] += 1

        #############################################################################
        
        outname = self.folder + self.sample + "_expected_variant_coverage_per_gene.tsv"
        
        with open(outname, "w") as out:
            
            header = ["gene", "discoverable", "total", "ratio"]
            header = "\t".join([str(x) for x in header])
            out.write(header + "\n")
            
            for gene in self.all_genes:
                total = self.gene_variants[gene]["total"]
                discoverable = self.gene_variants[gene]["discoverable"]
                
                try:
                    disco_ratio = round((discoverable/float(total)),2) *100
                except:
                    disco_ratio = 0
                
                if total > 0:
                    result = [gene, discoverable, total, disco_ratio]
                    result = "\t".join([str(x) for x in result])
                    out.write(result + "\n")
            
