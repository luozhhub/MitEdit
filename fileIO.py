import os
import sys
import pandas as pd
import pysam

class bamIO():
    def __init__(self, bam_file=None):
        self.bam = bam_file
        #read bam file
        self.samfile = pysam.AlignmentFile(self.bam, "rb")
        #check index file
        if not self.samfile.check_index():
            print("bam index file missing!")
            exit(1)
        
    def getBasicAlignInfo(self, chr="chrM"):
        """
        return read count mapped and unmapped. and information in each chromasome
        each chr information is like: IndexStats(contig='chr1', mapped=67654, unmapped=255, total=67909)
        
        contigs are sorted in decoy
        output is a csv file "mapping_stats.csv"
        """
        mapped_reads = self.samfile.mapped
        ummapped_reads = self.samfile.unmapped
        
        chr_mapped = self.samfile.get_index_statistics()
        list_chr = []
        decoy = {'contig': 'decoy', 'mapped': 0, 'unmapped': 0, 'total': 0}
        #convert to dict and add information to decoy
        for one_chr in chr_mapped:
            #chro, mapped, unmmaped = one_chr[0], one_chr[1], one_chr[2]
            one_d = dict(one_chr._asdict())
            if not one_d["contig"].startswith("chr"):
                decoy["mapped"] += one_d["mapped"]
                decoy["unmapped"] += one_d["unmapped"]
                decoy["total"] += one_d["total"]
                continue
            list_chr.append(one_d)
        list_chr.append(decoy)
        #construct dataframe
        df = pd.DataFrame(list_chr)
        states_file = os.path.join(os.path.dirname(self.bam), "mapping_stats.csv")
        df.to_csv(states_file, index=False)
        print("finish basic information statistic!")
        
    def getMTReads(self):
        """
        extract nitochondria read, properly paired read
        """
        #get MT read
        MTSam = self.samfile.fetch('chrM')
        #print(MTSam.get_reference_name)
        #output file
        MTSam_file = os.path.join(os.path.dirname(self.bam), "MT.bam")        
        pairedreads = pysam.AlignmentFile(MTSam_file, "wb", template=self.samfile)
        
        #get reference id
        tid = self.samfile.get_tid("chrM")
        #filter, properly paired, 
        for read in MTSam:
            chr = read.reference_id
            rnext = read.next_reference_id
            #condition is properly paied, chrom == "chrM", rnext=="chrM"
            if read.is_proper_pair and chr==tid and rnext==tid: 
                pairedreads.write(read)
        return(MTSam_file)
        
        
    def markduplicate(self):
        pass
        
        
        
if __name__ == "__main__":
    bamio = bamIO(bam_file = "/home/zyang/Project/mitochondria/pnas_data/luo_pipeline/work_test/bam/ERR452358/ERR452358.sort.bam")
    #bamio.getBasicAlignInfo()
    bamio.getMTReads()    
