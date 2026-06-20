
from sys import argv
import os
import multiprocessing

raw_bam,ref_file,threads=argv[1:]


os.system("blasr --nproc 60 --bam --out %s.blasr.bam %s %s"%(raw_bam[:-4],raw_bam,ref_file))


os.system("samtools view -@ 50 %s.blasr.bam |awk '$2 == 16 {print $1> \"%s.blasr.16.txt\"} $2 == 0 {print $1> \"%s.blasr.0.txt\"}'"%(raw_bam[:-4],raw_bam[:-13],raw_bam[:-13]))

print("zmw(0/16) number done!")


os.system("samtools view -H %s > %s.p.sam"%(raw_bam,raw_bam[:-4]))
os.system("samtools view -H %s > %s.n.sam"%(raw_bam,raw_bam[:-4]))


os.system("time samtools view %s | awk 'NR==FNR{a[$1]; next} $1 in a' %s.blasr.16.txt - >> %s.n.sam"%(raw_bam,raw_bam[:-13],raw_bam[:-4]))
os.system("time samtools view %s | awk 'NR==FNR{a[$1]; next} $1 in a' %s.blasr.0.txt - >> %s.p.sam"%(raw_bam,raw_bam[:-13],raw_bam[:-4]))

print("postive and negative strain done!")


#CCS
def run_ccs(name):
    os.system("samtools view -bS %s -o %s.bam"%(name,name[:-4]))
    #plz pay attention to the number of thread
    os.system("ccs --min-length 2800 --max-length 4000 -j %s --min-passes 3 %s.bam %s_ccs_pass3.bam"%(threads,name[:-4],name[:-4]))
    os.system("samtools view -h -@ 50 -o %s_ccs_pass3.sam %s_ccs_pass3.bam"%(name[:-4],name[:-4]))


p_n=[raw_bam[:-4]+".p.sam",raw_bam[:-4]+".n.sam"]
multiple_p=multiprocessing.Pool(2)
for file in p_n:
    ress=multiple_p.apply_async(run_ccs,args=(file,))
     
print('Waiting for all subprocesses done...')
multiple_p.close()
multiple_p.join()
print('All subprocesses done.')