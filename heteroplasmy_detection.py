import sys,os
import pandas as pd

from create_sample import create_sample_list

from subprocess import *

class Basic(object):
    def __init__(self):
        pass

    @staticmethod
    def run(cmd=None, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir, stdout=PIPE)
        p.wait()
        return p

class heteroplasmy(object):
    def __init__(self, work_dir="", sample_list_file=None, threads=None):
        """
        1. work_dir: the path should be a new dir, contain all result
        the structure of work dir:
            ---work dir
           |__bam
           |__mpile
           |__pbs
        
        2. sample_list_file: this file should contain samples by rows
        the sample_list_file format:
            sample_id\tfastq_1\tfastq_2
            
        3. threads is the cpu core used in bwa
        """
        if work_dir == "":
            print("work path should be setting!")
            exit(1)
        #parameter    
        self.prefix= work_dir
        if not os.path.exists(self.prefix):
            os.makedirs(self.prefix)
        self.sample_list_file = sample_list_file
        self.threads= threads        
        
        #tools
        self.bwa = "/home/nazhang/luozhihui/software/bwa/bwa"
        self.samtools = "/home/zxchen/anaconda3/bin/samtools"
        self.bbmap = "/home/zhluo/miniconda3/bin/filterbyname.sh"
        self.het_raw="/home/zyang/software/heteroplasmy.pyflow/heteroplasmy.mle.py"
        self.het_filter="/home/zyang/software/heteroplasmy.pyflow/heteroplasmy.mle.filter.py"
        #human reference hg19
        self.ref_genome = "/home/zyang/mit/ref/chrM.fa"
        self.whole_genome = "/home/zhluo/Project/mitch/autism/ref/GRCh37.p13.genome.fa"
        
        #output bam dir
        self.bam_dir = os.path.join(self.prefix, "bam")
        if not os.path.exists(self.bam_dir):
            os.makedirs(self.bam_dir)
        #output mpileup file dir
        self.mpile_dir = os.path.join(self.prefix, "mpile") 
        if not os.path.exists(self.mpile_dir):
            os.makedirs(self.mpile_dir)
        #pbs dir
        self.pbs_dir = os.path.join(self.prefix, "pbs") 
        if not os.path.exists(self.pbs_dir):
            os.makedirs(self.pbs_dir)
            
    def read_sample_file(self):
        """
        every sample should be put in one line,
        format:
        sample_id\tfastq_1\tfastq_2 
        """
        df = pd.read_csv(self.sample_list_file, sep="\t", names=["sampleID", "fastq_1", "fastq_2"])
        return(df)
        
    def process(self, sampleID=None, fastq_1=None, fastq_2=None):
        oneSample_dir = os.path.join(self.bam_dir, sampleID)
        if not os.path.exists(oneSample_dir):
            os.makedirs(oneSample_dir)
        log_file = os.path.join(oneSample_dir, sampleID + ".log.txt")
        #step1. bwa mapping to whole genome
        sort_bam = os.path.join(oneSample_dir, sampleID + ".sort.bam")
        cmd = """%s mem -R '@RG\\tID:GRCH38\\tSM:%s\\tLB:\\tPL:ILLUMINA' -t %s %s %s %s|%s view -@ %s -Shu -|%s sort -@ %s -o %s - > %s;\n""" % \
                (self.bwa, sampleID, self.threads, self.whole_genome, fastq_1, fastq_2, self.samtools, self.threads,self.samtools, self.threads,\
                 sort_bam, log_file)
        ##step1.1 bwa index for bam file
        cmd += """%s index %s;\n""" % (self.samtools, sort_bam)
        
        #step2. filter out read mapped in bam file
        MT_bam = os.path.join(oneSample_dir, sampleID + ".MT.bam")
        cmd += "%s view -bh %s chrM -o %s;\n" %(self.samtools, sort_bam, MT_bam)

        #step3. obtain unique mapped reads
        MT_bam_unique = os.path.join(oneSample_dir, sampleID + ".unique.MT.bam")
        MT_bam_noneUnique = os.path.join(oneSample_dir, sampleID + ".noneUnique.MT.bam")
        cmd += "%s view -h %s | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > %s;" % \
                (self.samtools, MT_bam, MT_bam_unique)
        cmd += "%s view -h %s | grep -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > %s;" % \
                (self.samtools, MT_bam, MT_bam_noneUnique)
        
        #step4. mpileup 
        p_sort_bam = os.path.join(oneSample_dir, sampleID + ".p.sort.bam")
        cmd += "%s sort -@ 7 -o %s %s;\n" % (self.samtools, p_sort_bam, MT_bam_unique)
        mpile_file = os.path.join(self.mpile_dir, sampleID + ".mpile.file")
        cmd += "%s mpileup -B -Q 30 -d 1000000 -L 10000 -f %s %s > %s;" %(self.samtools, self.ref_genome, p_sort_bam, mpile_file)

        #step5. call heteroplasmy
        # heteroplasmy
        heteroplasmy_raw = os.path.join(oneSample_dir, sampleID + ".mp.raw")
        cmd += '%s -i %s -o %s;\n' %(self.het_raw, mpile_file, heteroplasmy_raw)
        heteropasmy = os.path.join(self.mpile_dir, sampleID + ".heteroplasmy.txt")
        cmd += 'python2 %s --loose %s --chi %s -d %s --mle %s -i %s -o %s' \
                             %(self.het_filter, "0.003,0.003", "0.0", "1000", "0", heteroplasmy_raw, heteropasmy)
        
        return(cmd)
        
    def alignmrnt_parse(self, alignment_file):
        """
        alignment_file is sort bam file
        """
        pass

    def main(self):
        df_sample = self.read_sample_file()
        for idx,row in df_sample.iterrows():
            sampleID = row["sampleID"]
            fastq_1 = row["fastq_1"]
            fastq_2 = row["fastq_2"]
            
            cmd = self.process(sampleID=sampleID, fastq_1=fastq_1, fastq_2=fastq_2)
            handle = open(os.path.join(self.pbs_dir , sampleID + ".pbs"), "w")
            handle.write(cmd)
            handle.close()
        print("finish pbs job!")
        
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
            
            
        
    
if __name__ == "__main__":
    """
    todo, I may have to do quality trim for fastq file tommorrow
    """
    #create work dir
    work_dir = "/home/zyang/Project/mitochondria/pnas_data/luo_pipeline/work_test"
    if not os.path.exists(work_dir):
            os.makedirs(work_dir)
    #step1. create fastq input file
    fastq_dir = "/home/zyang/Project/mitochondria/pnas_data/fastq"
    create_sample_file = create_sample_list(sample_dir=fastq_dir, work_dir=work_dir)
    sample_list_file = create_sample_file.create_sample_file()
    #step2. create pbs
    call_mut = heteroplasmy(work_dir=work_dir, sample_list_file=sample_list_file, threads=7)
    call_mut.main()
    
    #step 3: submit qsub
    call_mut.submit(number="last_one", mem="15G")