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
        self.sample_meta_infor = {}
        
    def getBasicAlignInfo(self, chr="chrM"):
        """
        return read count mapped and unmapped. and information in each chromasome
        each chr information is like: IndexStats(contig='chr1', mapped=67654, unmapped=255, total=67909)
        
        contigs are sorted in decoy
        output is a csv file "mapping_stats.csv"
        """
        mapped_reads = self.samfile.mapped
        ummapped_reads = self.samfile.unmapped
        self.sample_meta_infor["mapped_reads"] = mapped_reads
        
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
        extract nitochondria read
        """
        #get MT read
        
        MTSam = self.samfile.fetch('chrM')
        return(MTSam)
        
        
    def properly_paired(self, samfile=None):
        """
        get properly paired read rnext "="
        """
        #output file
        MTSam_file = os.path.join(os.path.dirname(self.bam), "MT.bam")        
        pairedreads = pysam.AlignmentFile(MTSam_file, "wb", template=self.samfile)
        
        #get reference id
        tid = self.samfile.get_tid("chrM")
        #filter, properly paired,
        MTSam = samfile 
        properly_n, diff_chr_n = 0, 0
        MT_n = 0
        for read in MTSam:
            MT_n += 1
            chr = read.reference_id
            rnext = read.next_reference_id
            #condition is properly paied, chrom == "chrM", rnext=="chrM"
            if read.is_proper_pair:
                properly_n = properly_n + 1
                if chr==tid and rnext==tid:
                    diff_chr_n = diff_chr_n + 1 
                    pairedreads.write(read)
        self.sample_meta_infor["MT_reads"] = MT_n
        self.sample_meta_infor["properly_paired"] = properly_n
        self.sample_meta_infor["sameChrom_paired"] = diff_chr_n
        return(MTSam_file)
        
    def mapped_read(self, input_bam=None):
        samfile = pysam.AlignmentFile(input_bam, "rb")
        mapped_reads = samfile.mapped
        return(mapped_reads)
        
        
    def limit_mismatch(self, input_bam=None, output_bam=None, mismatch_num=2):
        """
        limit mismatch number of each read.
        default is 2
        only mapped read has NM field. so this step is after MT bam extracting
        """
        bamfile = pysam.AlignmentFile(input_bam, "rb")
        if not bamfile.check_index():
            print("bam index file missing!")
            exit(1)
        filter_bam = output_bam
        pairedreads = pysam.AlignmentFile(filter_bam, "wb", template=bamfile)
        before, after = 0, 0
        for read in bamfile:
            tag = read.get_tag("NM", with_value_type=False)
            before += 1
            if (tag < mismatch_num):
                pairedreads.write(read)
                after += 1
        pairedreads.close()
        bamfile.close()
        self.sample_meta_infor["befor_mismatch_filter"] = before
        self.sample_meta_infor["after_mismatch_filter"] = after
        
              
        
    def run(self):
        """
        all bam process step
        """
        #step1. get MT read, properly paired read, fq1 and fq2 both mapped to chrM read
        self.getBasicAlignInfo()
        MTSam = self.getMTReads()
        filterSam = self.properly_paired(samfile=MTSam)
        print(self.sample_meta_infor)
        return(filterSam)
        
    
        
if __name__ == "__main__":
    bamio = bamIO(bam_file = "/home/zyang/Project/mitochondria/pnas_data/luo_pipeline/work_test/bam/ERR452358/ERR452358.sort.bam")    
    #bamio.run()
    bamio.limit_mismatch(input_bam="/home/zyang/Project/mitochondria/pnas_data/luo_pipeline/work_test/bam/ERR452358/ERR452358.localRealign.bam")
        
