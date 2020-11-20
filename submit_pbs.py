import sys
import os
import pandas as pd
import argparse

from basic import Basic
from create_sample import create_sample_list

class job_submit():
    def __init__(self, sample_list_file=None, threads=None, pbs_dir= None):
        self.sample_list_file = sample_list_file
        self.threads = threads
        self.pbs_dir = pbs_dir

    def submit(self, number=None, mem="15G"):
        """
        number: how many jobs to submit:
        "all", "last_one"
        """
        handle = open(self.sample_list_file, "r")
        samples = handle.readlines()
        if number == "last_one":
            last_one = samples[-1].split("\t")[0]
            last_one_pbs = os.path.join(self.pbs_dir, last_one + ".pbs")
            cmd = "qsub -l nodes=1:ppn=%s -l mem=%s %s" %(self.threads, mem, last_one_pbs)
            p = Basic.run(cmd, wkdir=self.pbs_dir)
            print(p.stdout.read())
            
        if number == "all":
            for one_sample in samples:
                one = one_sample[-1].split("\t")[0]
                one_pbs = os.path.join(self.pbs_dir, one + ".pbs")
                cmd = "qsub -l nodes=1:ppn=%s -l mem=%s %s" %(self.threads, mem, one_pbs)
                p = Basic.run(cmd, wkdir=self.pbs_dir)
                print(p.stdout.read())
        print("finish submit pbs!")

class create_pbs():
    def __init__(self, sample_list_file=None, threads=None, pbs_dir= None, config_file=None, work_dir = None, fastq_dir = None):
        self.sample_list_file = sample_list_file
        self.pbs_dir = pbs_dir
        self.config_file = config_file
        self.work_dir = work_dir
        #self.fastq_dir = fastq_dir
    def read_sample_file(self):
        """
        every sample should be put in one line,
        format:
        "sampleID\tfastq_1\tfastq_2\n" 
        """
        df = pd.read_csv(self.sample_list_file, sep="\t", names=["sampleID", "fastq_1", "fastq_2"])
        return(df)
        
    def create_pbs(self):
        df_sample = self.read_sample_file()
        for idx,row in df_sample.iterrows():
            sampleID = row["sampleID"]
            fastq_1 = row["fastq_1"]
            fastq_2 = row["fastq_2"]
            
            #cmd = self.process(sampleID=sampleID, fastq_1=fastq_1, fastq_2=fastq_2)
            cmd = "python /home/zyang/software/MitEdit/heteroplasmy_detection.py --sampleID %s --fastq_1 %s --fastq_2 %s --config %s --work_dir %s"% (sampleID, fastq_1, fastq_2, self.config_file, self.work_dir)
            handle = open(os.path.join(self.pbs_dir , sampleID + ".pbs"), "w")
            handle.write(cmd)
            handle.close()
        print("create all pbs!")
        
def run(prog=None, args=None):
    usage = "usage: python submit_pbs.py \n\
            submit job to call heteroplasmy for each sample\n\
            Author: Zhihui Luo \n\
            last update 11/19 2020"
    
    parser = argparse.ArgumentParser(description = usage)
    parser.add_argument("--config",type=str,required=False,default='/home/zyang/software/MitEdit/config.txt', help='config file containing setting in MitEdit')
    parser.add_argument("--work_dir",type=str,required=False,default="/home/zyang/Project/mitochondria/pnas_data/luo_pipeline/work_test", help='all output file dir')
    parser.add_argument("--fastq_dir",type=str,required=False,default="/home/zyang/Project/mitochondria/pnas_data/fastq", help='dir containing the input fastq')
    opts = parser.parse_args(args)
        
    #argument = Basic.read_arguments(opts.config)
    
    #work directory
    work_dir = opts.work_dir
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
        
    #fastq dir
    fastq_dir = opts.fastq_dir
    
    #pbs_dir
    pbs_dir = os.path.join(work_dir, "pbs") 
    if not os.path.exists(pbs_dir):
        os.makedirs(pbs_dir)
    
    #step1. create fastq input file
    create_sample_file = create_sample_list(sample_dir=fastq_dir, work_dir=work_dir)
    sample_list_file = create_sample_file.create_sample_file()
    
    #step2. create pbs
    pbs = create_pbs(sample_list_file=sample_list_file, threads=None, pbs_dir= pbs_dir, config_file=opts.config, work_dir = opts.work_dir)
    pbs.create_pbs()
    
    #step3. submit pbs
    #job = job_submit(sample_list_file=sample_list_file, threads=None, pbs_dir= None)        
    #job.submit(number="last_one", mem="15G")     
        
if __name__ == "__main__":
    run(prog=sys.argv[0:1], args=sys.argv[1:])