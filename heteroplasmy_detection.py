import sys,os
import pandas as pd
import argparse

from create_sample import create_sample_list
from fileIO import bamIO
from basic import Basic



class heteroplasmy(object):
    def __init__(self, work_dir=None, threads=None):
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
        
        self.threads= threads        
        
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
        #tmp dir    
        self.tmp_dir = os.path.join(self.prefix, "tmp") 
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
            
    
        
    def process(self, sampleID=None, fastq_1=None, fastq_2=None):
        oneSample_dir = os.path.join(self.bam_dir, sampleID)
        if not os.path.exists(oneSample_dir):
            os.makedirs(oneSample_dir)
        log_file = os.path.join(oneSample_dir, sampleID + ".log.txt")
        #step1. bwa mapping to whole genome
        sort_bam = os.path.join(oneSample_dir, sampleID + ".sort.bam")
        cmd = """%s mem -R '@RG\\tID:GRCH37\\tSM:%s\\tLB:\\tPL:ILLUMINA' -t %s %s %s %s|%s view -@ %s -Shu -|%s sort -@ %s -o %s - > %s;\n""" % \
                (self.bwa, sampleID, self.threads, self.whole_genome, fastq_1, fastq_2, self.samtools, self.threads,self.samtools, self.threads,\
                 sort_bam, log_file)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        ##step1.1 bwa index for bam file
        cmd = """%s index %s;\n""" % (self.samtools, sort_bam)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        
        #step2. filter out read mapped in bam file
        MT_bam = os.path.join(oneSample_dir, sampleID + ".MT.bam")
        cmd = "%s view -bh %s chrM -o %s;\n" %(self.samtools, sort_bam, MT_bam)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        
        #step3. obtain unique mapped reads
        MT_bam_unique = os.path.join(oneSample_dir, sampleID + ".unique.MT.bam")
        MT_bam_noneUnique = os.path.join(oneSample_dir, sampleID + ".noneUnique.MT.bam")
        cmd = "%s view -h %s | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > %s;" % \
                (self.samtools, MT_bam, MT_bam_unique)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        cmd = "%s view -h %s | grep -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > %s;" % \
                (self.samtools, MT_bam, MT_bam_noneUnique)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        #step4. mpileup 
        p_sort_bam = os.path.join(oneSample_dir, sampleID + ".p.sort.bam")
        cmd = "%s sort -@ 7 -o %s %s;\n" % (self.samtools, p_sort_bam, MT_bam_unique)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        mpile_file = os.path.join(self.mpile_dir, sampleID + ".mpile.file")
        cmd = "%s mpileup -B -Q 30 -d 1000000 -L 10000 -f %s %s > %s;" %(self.samtools, self.ref_genome, p_sort_bam, mpile_file)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        
        #step5. call heteroplasmy
        # heteroplasmy
        heteroplasmy_raw = os.path.join(oneSample_dir, sampleID + ".mp.raw")
        cmd = '%s -i %s -o %s;\n' %(self.het_raw, mpile_file, heteroplasmy_raw)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        heteropasmy = os.path.join(self.mpile_dir, sampleID + ".heteroplasmy.txt")
        cmd = 'python2 %s --loose %s --chi %s -d %s --mle %s -i %s -o %s' \
                             %(self.het_filter, "0.003,0.003", "0.0", "1000", "0", heteroplasmy_raw, heteropasmy)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        print("finish detetcting!")
        
    def new_process(self, sampleID=None, fastq_1=None, fastq_2=None):
        oneSample_dir = os.path.join(self.bam_dir, sampleID)
        if not os.path.exists(oneSample_dir):
            os.makedirs(oneSample_dir)
        log_file = os.path.join(oneSample_dir, sampleID + ".log.txt")
        #step1. bwa mapping to whole genome
        sort_bam = os.path.join(oneSample_dir, sampleID + ".sort.bam")
        cmd = """%s mem -R '@RG\\tID:GRCH37\\tSM:%s\\tLB:mitochondira\\tPL:ILLUMINA' -t %s %s %s %s|%s view -@ %s -Shu -|%s sort -@ %s -o %s - > %s;\n""" % \
                (self.bwa, sampleID, self.threads, self.whole_genome, fastq_1, fastq_2, self.samtools, self.threads,self.samtools, self.threads,\
                 sort_bam, log_file)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        cmd = """%s index %s;\n""" % (self.samtools, sort_bam)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        
        #step2. filter out read mapped in bam file
        bamio = bamIO(bam_file = "/home/zyang/Project/mitochondria/pnas_data/luo_pipeline/work_test/bam/ERR452358/ERR452358.sort.bam")
        MTbam = bamio.run() 
        cmd = """%s index %s;\n""" % (self.samtools, MTbam)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        
        
        #step3. mark duplicate
        marked_duplicated_bam = os.path.join(oneSample_dir, sampleID + ".mkdup.bam")
        mkmatrix = os.path.join(oneSample_dir, sampleID + ".mkdup.matrix.txt") 
        cmd = self.markduplicate(input_bam = MTbam, marked_duplicated_bam=marked_duplicated_bam, marked_dup_metrics_txt=mkmatrix)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        bamio.sample_meta_infor["befor_mkduped_read"] = bamio.mapped_read(input_bam = MTbam)
        cmd = """%s index %s;\n""" % (self.samtools, marked_duplicated_bam)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        bamio.sample_meta_infor["after_mkduped_read"] = bamio.mapped_read(input_bam = marked_duplicated_bam)
        
        
        #step4. local realignment
        localRealign_bam = os.path.join(oneSample_dir, sampleID + ".localRealign.bam")
        cmd = self.localRealignment(input_bam=marked_duplicated_bam, output_bam=localRealign_bam)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        cmd = """%s index %s;\n""" % (self.samtools, localRealign_bam)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        bamio.sample_meta_infor["after_local_realign"] = bamio.mapped_read(input_bam = localRealign_bam)
        
        #step5. mismatch filter
        mismatch_bam = os.path.join(oneSample_dir, sampleID + ".mismatchFilter.bam")
        bamio.limit_mismatch(input_bam=localRealign_bam, output_bam=mismatch_bam, mismatch_num=1)
        cmd = """%s index %s;\n""" % (self.samtools, localRealign_bam)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        
        
        #step5. mpileup 
        p_sort_bam = os.path.join(oneSample_dir, sampleID + ".p.sort.bam")
        MT_bam_unique = os.path.join(oneSample_dir, sampleID + ".unique.MT.bam")
        cmd = "%s view -h %s | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > %s;" % \
                (self.samtools, mismatch_bam, MT_bam_unique)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        cmd = "%s sort -@ 7 -o %s %s;\n" % (self.samtools, p_sort_bam, MT_bam_unique)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        mpile_file = os.path.join(self.mpile_dir, sampleID + ".mpile.file")
        cmd = "%s mpileup -B -Q 20 -q 20 -d 1000000 -L 10000 -f %s %s > %s;" %(self.samtools, self.ref_genome, MT_bam_unique, mpile_file)
        Basic.run(cmd=cmd, wkdir=self.prefix)

        #step5. call heteroplasmy
        # heteroplasmy
        heteroplasmy_raw = os.path.join(oneSample_dir, sampleID + ".mp.raw")
        cmd = '%s -i %s -o %s;\n' %(self.het_raw, mpile_file, heteroplasmy_raw)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        heteropasmy = os.path.join(self.mpile_dir, sampleID + ".heteroplasmy.txt")
        cmd = 'python2 %s --loose %s --chi %s -d %s --mle %s -i %s -o %s' \
                             %(self.het_filter, "0.003,0.003", "0.0", "1000", "0", heteroplasmy_raw, heteropasmy)
        Basic.run(cmd=cmd, wkdir=self.prefix)
        print("finish heteroplasmy detecting!")
        print(bamio.sample_meta_infor)
        
    def argument_parse(self, config_file=None):
        """
        config_file: including all tools path
        this function set attribute to self, all tools in config file will be set
        #tools
        bwa = /home/zhluo/miniconda3/bin/bwa
        samtools = /home/zxchen/anaconda3/bin/samtools
        picard = /home/zhluo/Software/picard.jar
        bbmap = /home/zhluo/miniconda3/bin/filterbyname.sh
        het_raw=/home/zyang/software/heteroplasmy.pyflow/heteroplasmy.mle.py
        het_filter=/home/zyang/software/heteroplasmy.pyflow/heteroplasmy.mle.filter.py
        bamleftalign= /home/zhluo/miniconda3/bin/bamleftalign
        #human reference hg19
        ref_genome = /home/zyang/mit/ref/chrM.fa
        whole_genome = /home/zhluo/Project/mitch/autism/ref/GRCh37.p13.genome.fa
        
        """
        arguments = Basic.read_arguments(arg_file=config_file)
        for key in arguments:
            setattr(self, key, arguments[key])        
        
    def markduplicate(self, input_bam=None, marked_duplicated_bam=None, marked_dup_metrics_txt=None):
        """
        picard markduplicate
        """
        cmd="java -Djava.io.tmpdir=%s -jar -Xms1024m -Xmx10240m %s MarkDuplicates I=%s O=%s M=%s VALIDATION_STRINGENCY=LENIENT ASSUME_SORT_ORDER=queryname REMOVE_DUPLICATES=true;\n"\
            %(self.tmp_dir, self.picard, input_bam, marked_duplicated_bam, marked_dup_metrics_txt)
        return(cmd)
        
    def localRealignment(self, input_bam=None, output_bam=None):
        """
        this function is a wrap for bamleftalign.
        local realignment and recalibrate base quality
        """
        cmd = "%s -c -f %s < %s | %s calmd -EArb - %s >%s;\n" % (self.bamleftalign, self.ref_genome, input_bam,
        self.samtools,  self.ref_genome, output_bam)
        return(cmd)


        
def run(prog=None, args=None):
    usage = """usage: 
python /home/zyang/software/MitEdit/heteroplasmy_detection.py --sampleID [] --fastq_1 [] --fastq_2 [] --config [] --work_dir []
call heteroplasmy for each sample
Author: Zhihui Luo
last update 11/19 2020"""
    #parser = argparse.ArgumentParser(, description = usage)
    print(args)
    parser = argparse.ArgumentParser(prog = prog, usage=usage)
    parser.add_argument("--sampleID",type=str,required=True, help='sample id of fastq file')
    parser.add_argument("--fastq_1",type=str,required=True, help='input fastq 1')
    parser.add_argument("--fastq_2",type=str,required=True, help='input fastq 2')
    parser.add_argument("--config",type=str,required=True, help='config file containing all excute file path')
    parser.add_argument("--work_dir",type=str,required=True, help='output dir')
    opts = parser.parse_args(args)
        
    #initial heteroplasmy class
    hetero = heteroplasmy(work_dir=opts.work_dir, threads=7)
    hetero.argument_parse(config_file=opts.config)
    hetero.new_process(sampleID=opts.sampleID, fastq_1=opts.fastq_1, fastq_2=opts.fastq_2)


if __name__ == "__main__":
    print(sys.argv[0])
    run(prog=sys.argv[0:1], args=sys.argv[1:])        