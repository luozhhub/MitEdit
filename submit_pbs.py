import sys
import os

from basic import Basic

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
        
if __name__ == "__main__":
    job = job_submit(sample_list_file=None, threads=None, pbs_dir= None)
    call_mut.submit(number="last_one", mem="15G") 