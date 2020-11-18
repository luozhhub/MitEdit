import pandas as pd
import sys,os

class create_sample_list():
    def __init__(self, sample_dir=None, work_dir=None):
        self.sample_dir = sample_dir
        self.sample_file = os.path.join(work_dir, "all_sample.txt")
        
    def create_sample_file(self):
        """
        This function is suit for data downloaded from SRA
        
        file name example:
        ERR452355_1.fastq.gz
        ERR452355_2.fastq.gz
        
        """
        sample_list = []
        for one_file in os.listdir(self.sample_dir):
            if ("fastq" in one_file) or ("fq" in one_file):
                array = one_file.split("_")
                #the third file of the sample sqlt 3
                if len(array) ==1:
                    continue
                sampleID = array[0]
                if sampleID not in sample_list:
                    sample_list.append(sampleID)
                    
        handle = open(self.sample_file, "w")
        for sampleID in sorted(sample_list):
            fastq_1 = os.path.join(self.sample_dir, sampleID + "_1.fastq.gz")
            fastq_2 = os.path.join(self.sample_dir, sampleID + "_2.fastq.gz")
            #check file exists
            if not os.path.exists(fastq_1) or not os.path.exists(fastq_2):
                print("the fastq sample name is not correct!")
                exit(1)
            #output 
            handle.write("%s\t%s\t%s\n" % (sampleID, fastq_1, fastq_2))
        handle.close()
        return(self.sample_file)
            
                            
                
            
        