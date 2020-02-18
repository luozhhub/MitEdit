#!/usr/bin/python3
"""

"""
import sys,os
import pandas as pd

class heteroplasmy():
    def __init__(self, prefix="/home/zyang/Project/mitochondria/blood_tissue/"):
        self.prefix= prefix
        #tools
        self.bwa = "/home/nazhang/luozhihui/software/bwa/bwa"
        self.samtools = "/home/zxchen/anaconda3/bin/samtools"
        self.bbmap = "/home/zhluo/miniconda3/bin/filterbyname.sh"
        #human reference
        self.ref_genome = "/home/zyang/mit/ref/chrM.fa"
        self.whole_genome = "/home/zhluo/Project/mitch/autism/ref/GRCh37.p13.genome.fa"
        #data dir
        self.fastq_dir = os.path.join(self.prefix, "fastq") #"/home/zyang/Project/mitochondria/CFS_ME/fastq"
        self.bam_dir = os.path.join(self.prefix, "bam") #"/home/zyang/Project/mitochondria/CFS_ME/bam"
        if not os.path.exists(self.bam_dir):
            os.makedirs(self.bam_dir)
        self.mpile_dir = os.path.join(self.prefix, "mpile") #"/home/zyang/Project/mitochondria/CFS_ME/mpile"
        if not os.path.exists(self.mpile_dir):
            os.makedirs(self.bam_dir)
        self.pbs_dir = os.path.join(self.prefix, "pbs") #"/home/zyang/Project/mitochondria/CFS_ME/pbs"
        if not os.path.exists(self.pbs_dir):
            os.makedirs(self.bam_dir)

#pair end

    def create_pipeline(self, threads=20):
        sample_list = []
        for sample in os.listdir(self.fastq_dir):
            if ".fastq" in sample:
                sample_name = sample.split("_")[0]
                if sample_name in sample_list:
                    continue
                sort_bam = os.path.join(self.bam_dir, sample_name + ".sort.bam")
                mpile_file = os.path.join(self.mpile_dir, sample_name + ".mpile.file")
                handle = open(os.path.join(self.pbs_dir , sample_name + ".pbs"), "w")
                #1. bwa mapping
                fastq_1 = os.path.join(self.fastq_dir, sample_name + "_1.fastq.gz")
                fastq_2 = os.path.join(self.fastq_dir, sample_name + "_2.fastq.gz")
                cmd = """%s mem -R '@RG\\tID:GRCH38\\tSM:%s\\tLB:\\tPL:ILLUMINA' -t %s %s %s %s|%s view -@ %s -Shu -|%s sort -@ %s -o %s - > %s;\n""" % \
                (self.bwa, sample_name, threads, self.whole_genome, fastq_1, fastq_2, self.samtools, threads,self.samtools, threads,sort_bam, os.path.join(self.bam_dir, sample_name + ".log.txt"))
                handle.write(cmd)
                #2. samtools index
                cmd = """%s index %s;\n""" % (self.samtools, sort_bam)
                handle.write(cmd)
                #3.
                MT_bam = os.path.join(self.bam_dir, sample_name + ".MT.bam")
                cmd = "%s view -bh %s chrM -o %s;\n" %(self.samtools, sort_bam, MT_bam)
                handle.write(cmd)
                #handle.write(cmd)
                MT_name = os.path.join(self.bam_dir, sample_name + "MT.name")
                cmd = """%s view %s |awk -F"\\t" '{print $1}'|sort|uniq -c|awk '{print $2}' > %s;\n""" % (self.samtools, MT_bam, MT_name)
                handle.write(cmd)
                #4.
                MT_fastq = os.path.join(self.bam_dir, sample_name + ".MT.fq")
                cmd = "%s in=%s out=%s names=%s include=t ow=t;\n" % (self.bbmap, os.path.join(self.fastq_dir, sample), MT_fastq, MT_name)
                handle.write(cmd)
                #5. 
                p_bam = os.path.join(self.bam_dir, sample_name + ".p.bam")
                cmd = "%s view -b -F 2304 %s > %s;\n" % (self.samtools, MT_bam,  p_bam)
                handle.write(cmd)
                #6. 
                p_sort_bam = os.path.join(self.bam_dir, sample_name + ".p.sort.bam")
                cmd = "%s sort -@ 12 -o %s %s;\n" % (self.samtools, p_sort_bam, p_bam)
                handle.write(cmd)
                #7.
                mpile_file = os.path.join(self.mpile_dir, sample_name + ".mpile.file")
                cmd = "%s mpileup -B -Q 30 -d 1000000 -L 10000 -f %s %s > %s;" %(self.samtools, self.ref_genome, p_sort_bam, mpile_file)
                
                handle.write(cmd)
                handle.close()
                sample_list.append(sample)
                
    def transfer_files(self):
        dest_dir = "/public/home/zhluo/project/mitch/blood_tissue/"
        cmd = "scp -r %s zhluo@211.69.141.130:%s" %(self.mpile_dir, dest_dir)
        os.system(cmd)
    
    def transfer_heteroresult(self):
        dest_dir = self.prefix
        cmd = "scp -r zhluo@211.69.141.130:%s %s" %("/public/home/zhluo/project/mitch/blood_tissue/mpile/heteroplasmy_result", dest_dir)
        os.system(cmd)
        
    def Parse_result(self):
        data_info = pd.read_csv("/home/zyang/Project/mitochondria/CFS_ME/cfs_me_metadata.txt", sep="\t", header=0)
        data_patient = data_info[data_info["filetype"] == "Patient"]
        data_control = data_info[data_info["filetype"] == "Control"]
        data_patient_day1 = data_patient[data_patient["date"] == 1]
        pd1_sra = list(data_patient_day1["Run"])
        data_patient_day2 = data_patient[data_patient["date"] == 2]
        pd2_sra = list(data_patient_day2["Run"])
        data_patient_day3 = data_patient[data_patient["date"] == 3]
        pd3_sra = list(data_patient_day3["Run"])
        data_patient_day7 = data_patient[data_patient["date"] == 7]
        pd7_sra = list(data_patient_day7["Run"])
        data_control_day1 = data_control[data_control["date"] == 1]
        cd1_sra = list(data_control_day1["Run"])
        data_control_day2 = data_control[data_control["date"] == 2]
        cd2_sra = list(data_control_day2["Run"])
        data_control_day3 = data_control[data_control["date"] == 3]
        cd3_sra = list(data_control_day3["Run"])
        data_control_day7 = data_control[data_control["date"] == 7]
        cd7_sra = list(data_control_day7["Run"])
        return [pd1_sra,pd2_sra,pd3_sra,pd7_sra,cd1_sra,cd2_sra,cd3_sra,cd7_sra]
        
        
    def statistic(self):
        result_dir = "/home/zyang/Project/mitochondria/CFS_ME/heteroplasmy_result"
        group_list = self.Parse_result()
        
        group_num = 0
        for one_group in group_list:
            group_num +=1
            flag = True
            sample_number = 0
            for one_file in os.listdir(result_dir):
                if "file_hmarker.txt" in one_file:
                    ids = one_file.split(".")[0]
                    if ids in one_group:
                        sample_number += 1
                        if flag is True:
                            df_tmp = pd.read_csv(os.path.join(result_dir, one_file), header=0, sep="\t", quotechar="'")
                            df_tmp["sample_id"] = ids
                            flag = False
                            continue
                        df = pd.read_csv(os.path.join(result_dir, one_file), header=0, sep="\t", quotechar="'")
                        df["sample_id"] = ids
                        df_tmp = df_tmp.append(df)
            df_tmp.to_csv("total_result_group_" + str(group_num)+".txt", sep="\t", header=True)
            result = df_tmp.groupby(["loc", "ref", "AssU", "Ualt"], as_index=False).count()
            result = result.sort_values(["heteroplasmy"], ascending =False)
            result = result[["loc", "ref", "AssU", "Ualt", "heteroplasmy"]]
            print(sample_number)
            result.to_csv("heteroplasmy_cfsME_group_" + str(group_num)+".result.txt", sep="\t", header=True)
                    
    def one_group_parse(self):
        result_dir = "/home/zyang/Project/mitochondria/blood_tissue/heteroplasmy_result"
        
        flag = True
        sample_number = 0
        for one_file in os.listdir(result_dir):
            if "file_hmarker.txt" in one_file:
                ids = one_file.split(".")[0]
                
                sample_number += 1
                if flag is True:
                    df_tmp = pd.read_csv(os.path.join(result_dir, one_file), header=0, sep="\t", quotechar="'")
                    df_tmp["sample_id"] = ids
                    flag = False
                    continue
                df = pd.read_csv(os.path.join(result_dir, one_file), header=0, sep="\t", quotechar="'")
                df["sample_id"] = ids
                df_tmp = df_tmp.append(df)
        df_tmp.to_csv("total_result.txt", sep="\t", header=True)
        result = df_tmp.groupby(["loc", "ref", "AssU", "Ualt"], as_index=False).count()
        result = result.sort_values(["heteroplasmy"], ascending =False)
        result = result[["loc", "ref", "AssU", "Ualt", "heteroplasmy"]]
        print(sample_number)
        result.to_csv("heteroplasmy.result.txt", sep="\t", header=True)
                

if __name__ == "__main__":
    hetero = heteroplasmy(prefix="/home/zyang/Project/mitochondria/blood_tissue/")
    #hetero.create_pipeline()
    #hetero.transfer_files()
    #hetero.transfer_heteroresult()
    hetero.one_group_parse()
    #hetero.Parse_result()
    #hetero.statistic()
    """
    for i in `ls ./*.pbs` ; do qsub -l nodes=1:ppn=20 $i ;done
    """