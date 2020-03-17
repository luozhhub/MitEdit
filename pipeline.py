#!/usr/bin/python3
"""

"""
import sys,os
import pandas as pd

class heteroplasmy():
    def __init__(self, prefix="/home/zyang/Project/mitochondria/CFS_ME/"):
        self.prefix= prefix
        #tools
        self.bwa = "/home/nazhang/luozhihui/software/bwa/bwa"
        self.samtools = "/home/zxchen/anaconda3/bin/samtools"
        self.bbmap = "/home/zhluo/miniconda3/bin/filterbyname.sh"
        #human reference
        self.ref_genome = "/home/zyang/mit/ref/chrM.fa"
        self.whole_genome = "/home/zhluo/Project/mitch/autism/ref/GRCh37.p13.genome.fa"
        #data dir
        #self.fastq_dir = os.path.join(self.prefix, "fastq") #"/home/zyang/Project/mitochondria/CFS_ME/fastq"
        self.fastq_dir = "/home/zyang/Project/mitochondria/CFS_cell_type/Fabien2699_Tcell_2015_10_09_FCB" 
        self.bam_dir = os.path.join(self.prefix, "bam") #"/home/zyang/Project/mitochondria/CFS_ME/bam"
        if not os.path.exists(self.bam_dir):
            os.makedirs(self.bam_dir)
        self.mpile_dir = os.path.join(self.prefix, "mpile") #"/home/zyang/Project/mitochondria/CFS_ME/mpile"
        if not os.path.exists(self.mpile_dir):
            os.makedirs(self.mpile_dir)
        self.pbs_dir = os.path.join(self.prefix, "pbs") #"/home/zyang/Project/mitochondria/CFS_ME/pbs"
        if not os.path.exists(self.pbs_dir):
            os.makedirs(self.pbs_dir)

#pair end

    def create_pipeline(self, threads=20):
        sample_list = []
        for sample in os.listdir(self.fastq_dir):
            if ".fastq" in sample:
                sample_name = sample.split(".")[0]
                if sample_name in sample_list:
                    continue
                sort_bam = os.path.join(self.bam_dir, sample_name + ".sort.bam")
                mpile_file = os.path.join(self.mpile_dir, sample_name + ".mpile.file")
                handle = open(os.path.join(self.pbs_dir , sample_name + ".pbs"), "w")
                #1. bwa mapping
                fastq_1 = os.path.join(self.fastq_dir, sample_name + ".fastq.gz")
                #fastq_2 = os.path.join(self.fastq_dir, sample_name + "_2.fastq.gz")
                cmd = """%s mem -R '@RG\\tID:GRCH38\\tSM:%s\\tLB:\\tPL:ILLUMINA' -t %s %s %s|%s view -@ %s -Shu -|%s sort -@ %s -o %s - > %s;\n""" % \
                (self.bwa, sample_name, threads, self.whole_genome, fastq_1, self.samtools, threads,self.samtools, threads,sort_bam, os.path.join(self.bam_dir, sample_name + ".log.txt"))
                handle.write(cmd)
                #2. samtools index
                cmd = """%s index %s;\n""" % (self.samtools, sort_bam)
                handle.write(cmd)
                #3.
                MT_bam = os.path.join(self.bam_dir, sample_name + ".MT.bam")
                cmd = "%s view -bh %s chrM -o %s;\n" %(self.samtools, sort_bam, MT_bam)
                handle.write(cmd)
                
                #3.5 obtain unique mapping read
                MT_bam_unique = os.path.join(self.bam_dir, sample_name + "unique.MT.bam")
                cmd = "%s view -h %s | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > %s;" % \
                (self.samtools, MT_bam, MT_bam_unique)
                handle.write(cmd)
                #handle.write(cmd)
                MT_name = os.path.join(self.bam_dir, sample_name + "MT.name")
                cmd = """%s view %s |awk -F"\\t" '{print $1}'|sort|uniq -c|awk '{print $2}' > %s;\n""" % (self.samtools, MT_bam_unique, MT_name)
                handle.write(cmd)
                #4.
                MT_fastq = os.path.join(self.bam_dir, sample_name + ".MT.fq")
                cmd = "%s in=%s out=%s names=%s include=t ow=t;\n" % (self.bbmap, os.path.join(self.fastq_dir, sample), MT_fastq, MT_name)
                handle.write(cmd)
                #5. 
                p_bam = os.path.join(self.bam_dir, sample_name + ".p.bam")
                cmd = "%s view -b -F 2304 %s > %s;\n" % (self.samtools, MT_bam_unique,  p_bam)
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
        dest_dir = "/public/home/zhluo/project/mitch/cfs_cell/mpile_02/"
        cmd = "scp -r %s zhluo@211.69.141.130:%s" %(self.mpile_dir, dest_dir)
        os.system(cmd)
    
    def transfer_heteroresult(self):
        dest_dir = "/home/zyang/Project/mitochondria/CFS_cell_type/Fabein2699_2015_06_12/"
        cmd = "scp -r zhluo@211.69.141.130:%s %s" %("/public/home/zhluo/project/mitch/cfs_cell/mpile_02/mpile/heteroplasmy_result", dest_dir)
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
        result_dir = "/home/zyang/Project/mitochondria/CFS_ME/pipeline_rmunique/heteroplasmy_result"
        group_list = self.Parse_result()
        
        group_num = 0
        groups = ["patient_day_1", "patient_day_2", "patient_day_3", "patient_day_7", "control_day_1", "control_day_2", "control_day_3", "control_day_7"]
        coverage = open("./depth_each_file.txt", "w")
        for one_group in group_list:
            
            flag = True
            sample_number = 0
            for one_file in os.listdir(result_dir):
                if "file_hmarker.txt" in one_file:
                    ids = one_file.split(".")[0]
                    if ids in one_group:
                        sample_number += 1
                        if flag is True:
                            df_tmp = pd.read_csv(os.path.join(result_dir, one_file), header=0, sep="\t", quotechar="'")
                            #df_tmp["sample_id"] = ids
                            df_tmp.insert(0, "sample_id", ids)
                            flag = False
                            continue
                        df = pd.read_csv(os.path.join(result_dir, one_file), header=0, sep="\t", quotechar="'")
                        #df["sample_id"] = ids
                        df.insert(0, "sample_id", ids)
                        df_tmp = df_tmp.append(df)
            df_tmp = df_tmp.sort_values(by=['loc'])
            sampe_name_group = groups[group_num]
            #print(df_tmp.loc[:, ["UA", "UT", "UC", "UG", "LA", "LT", "LC", "LG"]].sum(axis=1))
            #exit(1)
            df_tmp.insert(3, "total_read", list(df_tmp.loc[:, ["UA", "UT", "UC", "UG", "LA", "LT", "LC", "LG"]].sum(axis=1)))
            df_tmp.insert(4, "A_read", list(df_tmp.loc[:, ["UA","LA"]].sum(axis=1)))
            df_tmp.insert(5, "T_read", list(df_tmp.loc[:, ["UT","LT"]].sum(axis=1)))
            df_tmp.insert(6, "C_read", list(df_tmp.loc[:, ["UC","LC"]].sum(axis=1)))
            df_tmp.insert(7, "G_read", list(df_tmp.loc[:, ["UG","LG"]].sum(axis=1)))
            major_allele, minor_allele, major_frequence, minor_frequence = [], [], [], []
            for idx, row in df_tmp.iterrows():
                total = row["total_read"]
                A = row["A_read"]
                T = row["T_read"]
                C = row["C_read"]
                G = row["G_read"]
                arr = pd.Series({"A_read": A, "T_read":T, "C_read": C, "G_read": G})
                arr_nu = arr.sort_values(ascending=False).index
                major_allele.append(arr_nu[0][0])
                major_frequence.append(round(float(row[arr_nu[0]])/ total, 3))
                minor_allele.append(arr_nu[1][0])
                minor_frequence.append(round(float(row[arr_nu[1]])/ total, 3))
            df_tmp.insert(4, "major_allele", major_allele)
            df_tmp.insert(5, "major_frequence", major_frequence)
            df_tmp.insert(6, "minor_allele", minor_allele)
            df_tmp.insert(7, "minor_frequence", minor_frequence)  
            df_tmp = df_tmp[df_tmp["total_read"] >= 100]  
            df_tmp = df_tmp.sort_values(by=['minor_frequence', 'loc'], ascending=[False, True])
            df_tmp.to_csv("total_result_group_" + str(sampe_name_group)+".csv", sep=",", header=True)
            result = df_tmp.groupby(["loc", "ref", "AssU", "Ualt"], as_index=False).count()
            result = result.sort_values(["heteroplasmy"], ascending =False)
            result = result[["loc", "ref", "AssU", "Ualt", "heteroplasmy"]]
            print(sample_number)
            result.to_csv("heteroplasmy_cfsME_group_" + str(sampe_name_group)+".result.csv", sep=",", header=True)
            group_num +=1
            
            
            for one_file in os.listdir("/home/zyang/Project/mitochondria/CFS_ME/bam"):
                if "p.sort.bam" in one_file:
                    ids = one_file.split(".")[0]
                    if ids in one_group:
                        cmd = "samtools flagstat %s > %s" % (os.path.join("/home/zyang/Project/mitochondria/CFS_ME/bam", one_file), os.path.join("/home/zyang/Project/mitochondria/CFS_ME/bam", ids + ".flagstat.txt"))
                        os.system(cmd)
                        
                        out_file = open(os.path.join("/home/zyang/Project/mitochondria/CFS_ME/bam", ids + ".flagstat.txt"), "r")
                        lines = out_file.readlines()
                        for line in lines:
                            if "mapped (" in line:
                                array = line.split(" ")
                                coverage.write(sampe_name_group + "," +ids + "," + array[0] + "\n")
        coverage.close()
                        
    def cfs_cell(self):
        result_dir = "/home/zyang/Project/mitochondria/CFS_cell_type/Fabein2699_2015_06_25/heteroplasmy_result"
        start, end = 1, 99
        ctrl, patient = [], []
        for num in range(start, end + 1):
            if num % 2 != 0:
                ctrl.append(num)
            else:
                patient.append(num)
        sample_number = 0
        flag = True
        
        for one_file in os.listdir(result_dir):
            ids = one_file.split("_")[0].split("-")[1]
            if int(ids) in ctrl:
                if "file_hmarker.txt" in one_file:
                    sample_number += 1
                    if flag is True:
                        df_tmp = pd.read_csv(os.path.join(result_dir, one_file), header=0, sep="\t", quotechar="'")
                        #df_tmp["sample_id"] = ids
                        df_tmp.insert(0, "sample_id", ids)
                        flag = False
                        continue
                    df = pd.read_csv(os.path.join(result_dir, one_file), header=0, sep="\t", quotechar="'")
                    #df["sample_id"] = ids
                    df.insert(0, "sample_id", ids)
                    df_tmp = df_tmp.append(df)
                    #print (sample_number)
        df_tmp = df_tmp.sort_values(by=['loc'])
        #sampe_name_group = groups[group_num]
        #print(df_tmp.loc[:, ["UA", "UT", "UC", "UG", "LA", "LT", "LC", "LG"]].sum(axis=1))
        #exit(1)
        df_tmp.insert(3, "total_read", list(df_tmp.loc[:, ["UA", "UT", "UC", "UG", "LA", "LT", "LC", "LG"]].sum(axis=1)))
        df_tmp.insert(4, "A_read", list(df_tmp.loc[:, ["UA","LA"]].sum(axis=1)))
        df_tmp.insert(5, "T_read", list(df_tmp.loc[:, ["UT","LT"]].sum(axis=1)))
        df_tmp.insert(6, "C_read", list(df_tmp.loc[:, ["UC","LC"]].sum(axis=1)))
        df_tmp.insert(7, "G_read", list(df_tmp.loc[:, ["UG","LG"]].sum(axis=1)))
        major_allele, minor_allele, major_frequence, minor_frequence = [], [], [], []
        for idx, row in df_tmp.iterrows():
            total = row["total_read"]
            A = row["A_read"]
            T = row["T_read"]
            C = row["C_read"]
            G = row["G_read"]
            arr = pd.Series({"A_read": A, "T_read":T, "C_read": C, "G_read": G})
            arr_nu = arr.sort_values(ascending=False).index
            major_allele.append(arr_nu[0][0])
            major_frequence.append(round(float(row[arr_nu[0]])/ total, 3))
            minor_allele.append(arr_nu[1][0])
            minor_frequence.append(round(float(row[arr_nu[1]])/ total, 3))
        df_tmp.insert(4, "major_allele", major_allele)
        df_tmp.insert(5, "major_frequence", major_frequence)
        df_tmp.insert(6, "minor_allele", minor_allele)
        df_tmp.insert(7, "minor_frequence", minor_frequence)  
        df_tmp = df_tmp[df_tmp["total_read"] >= 100]  
        df_tmp = df_tmp.sort_values(by=['minor_frequence', 'loc'], ascending=[False, True])
        df_tmp.to_csv("result_ctrl" +".csv", sep=",", header=True)
        result = df_tmp.groupby(["loc", "ref", "AssU", "Ualt"], as_index=False).count()
        result = result.sort_values(["heteroplasmy"], ascending =False)
        result = result[["loc", "ref", "AssU", "Ualt", "heteroplasmy"]]
        print(sample_number)
        result.to_csv("heteroplasmy_cfs_cell_ctrl" +".result.csv", sep=",", header=True)
        #group_num +=1

    def R_script(self):
        handle = open("mito.pbs", "w")
        cmd = "/home/zhluo/miniconda3/bin/Rscript  /home/zyang/Project/mitochondria/CFS_cell_type/mitomut.R /home/zyang/Project/mitochondria/CFS_cell_type/Fabein2699_B_cell_2015_09_29/mpile/;"
        handle.write(cmd)
        cmd = "/home/zhluo/miniconda3/bin/Rscript  /home/zyang/Project/mitochondria/CFS_cell_type/mitomut.R /home/zyang/Project/mitochondria/CFS_cell_type/Fabien2699_NKCell_2015_10_05_FCA/mpile/;"
        handle.write(cmd)
        cmd = "/home/zhluo/miniconda3/bin/Rscript  /home/zyang/Project/mitochondria/CFS_cell_type/mitomut.R /home/zyang/Project/mitochondria/CFS_cell_type/Fabien2699_NKCell_2015_10_05_FCB/mpile/;"
        handle.write(cmd)
        cmd = "/home/zhluo/miniconda3/bin/Rscript  /home/zyang/Project/mitochondria/CFS_cell_type/mitomut.R /home/zyang/Project/mitochondria/CFS_cell_type/Fabien2699_Tcell_2015_10_09_FCA/mpile/;"
        handle.write(cmd)
        cmd = "/home/zhluo/miniconda3/bin/Rscript  /home/zyang/Project/mitochondria/CFS_cell_type/mitomut.R /home/zyang/Project/mitochondria/CFS_cell_type/Fabien2699_Tcell_2015_10_09_FCB/mpile/;"
        handle.write(cmd)
        handle.close()
                

if __name__ == "__main__":
    hetero = heteroplasmy(prefix="/home/zyang/Project/mitochondria/CFS_cell_type/Fabien2699_Tcell_2015_10_09_FCB")    
    #hetero.create_pipeline()
    #hetero.transfer_files()
    #hetero.transfer_heteroresult()
    #hetero.Parse_result()
    #hetero.statistic()
    #hetero.cfs_cell()
    """
    for i in "ls ./*.pbs" ; do qsub -l nodes=1:ppn=20 $i ;done
    """
    hetero.R_script() 
