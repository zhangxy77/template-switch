from sys import argv
import os
import re
import multiprocessing

raw_bam,ref_file,logfile,threads,=argv[1:]

os.system("blasr --nproc 30 --bam --out %s.blasr.bam %s %s"%(raw_bam[:-4],raw_bam,ref_file))

os.system("samtools sort -m 15G -n -o %s.blasr.sort.sam -O sam -@ 30 %s.blasr.bam"%(raw_bam[:-4],raw_bam[:-4]))
print("blasr.sort.sam is end")
os.system("samtools sort -m 15G -n -o %s.sort.sam -O sam -@ 30 %s"%(raw_bam[:-4],raw_bam))
print("sort.sam is end")

p_file=open("%s.p.sam"%(raw_bam[:-4]),"w")
n_file=open("%s.n.sam"%(raw_bam[:-4]),"w")
log_file=open(logfile,"w")

with open("%s.sort.sam"%(raw_bam[:-4]),"r") as sam_file:
    for hed in sam_file:
        if hed[0]!="m":
            p_file.write(hed)
            n_file.write(hed)
        else:
            break
def file(input_bam):
    with open(input_bam,"r") as bam:
        bam_dict={}
        count=0
        zmv1=""
        while True:
            try:
                l_1=next(bam).strip()
                if l_1[0]=="m":
                    hole_all=l_1.strip().split()[0]
                    hole_only=hole_all.split("/")[1]
                    if count>0 and hole_only!=zmv1:
                        yield bam_dict, zmv1
                        bam_dict={}
                        bam_dict[hole_all]=l_1
                        zmv1=hole_only
                    else:#
                        bam_dict[hole_all]=l_1
                        zmv1=hole_only
                        count+=1
            except StopIteration:
                yield bam_dict, zmv1
                print("Bam file is end")
                break
def b_file(input_blasr):
    with open(input_blasr,"r") as bam:
        bam_dict={}
        count=0
        zmv1=""
        while True:
            try:
                l_1=next(bam).strip()
                if l_1[0]=="m":
                    hole_all=l_1.split()[0]
                    hole_only=hole_all.split("/")[1]
                    if count>0 and hole_only!=zmv1:
                        yield bam_dict, zmv1
                        bam_dict={}
                        bam_dict[hole_all]=l_1.split()[1]
                        zmv1=hole_only
                    else:
                        bam_dict[hole_all]=l_1.split()[1]
                        zmv1=hole_only
                        count+=1
            except StopIteration:
                yield bam_dict, zmv1
                print("Blasr-bam file is end")
                break

mud_dict={"A":"T","T":"A","C":"G","G":"C"}
Reverse_complemrnt= lambda x:"".join([mud_dict[i] for i in x][::-1])


bam_file=file("%s.sort.sam"%(raw_bam[:-4]))
blasr_file=b_file("%s.blasr.sort.sam"%(raw_bam[:-4]))
c=0
while True:
    try:
        bam_subdict,zmv1=next(bam_file)
        blasr_subdict,zmv2=next(blasr_file)
        c+=1
        while zmv1!=zmv2:
            log_file.write(">"+zmv1+"\n")
            for ii in bam_subdict.keys():
                log_file.write(bam_subdict[ii]+"\n")
            bam_subdict,zmv1=next(bam_file)
        if zmv1==zmv2:
            for i in blasr_subdict:
                if blasr_subdict[i]=="16":
                    n_file.write(bam_subdict[i]+"\n")
                elif blasr_subdict[i]=="0":
                    p_file.write(bam_subdict[i]+"\n")
        else:
            print("error")
    except StopIteration:
        print("bam is end",zmv1,zmv2)
        break

p_file.close()
n_file.close()
log_file.close()

#for 1000bp:
def run_ccs(name):
    os.system("samtools view -bS %s.sam -o %s.bam"%(name[:-4],name[:-4]))
    #plz pay attention to the number of thread
    os.system("ccs --min-length 900 --max-length 2000 -j %s --min-passes 5  %s.bam %s_ccs.bam"%(threads,name[:-4],name[:-4]))
    os.system("samtools view -h -o %s_ccs.sam %s_ccs.bam"%(name[:-4],name[:-4]))

#for 500bp:
# def run_ccs(name):
#     os.system("samtools view -bS %s.sam -o %s.bam"%(name[:-4],name[:-4]))
#     #plz pay attention to the number of thread
#     os.system("ccs --min-length 400 --max-length 1000 -j %s --min-passes 5 %s.bam %s_ccs.bam"%(threads,name[:-4],name[:-4]))
#     os.system("samtools view -h -o %s_ccs.sam %s_ccs.bam"%(name[:-4],name[:-4]))

#for 250bp:
# def run_ccs(name):
#     os.system("samtools view -bS %s.sam -o %s.bam"%(name[:-4],name[:-4]))
#     #plz pay attention to the number of thread
#     os.system("ccs --min-length 200 --max-length 600 -j %s --min-passes 5  %s.bam %s_ccs.bam"%(threads,name[:-4],name[:-4]))
#     os.system("samtools view -h -o %s_ccs.sam %s_ccs.bam"%(name[:-4],name[:-4]))

p_n=[raw_bam[:-4]+".p.sam",raw_bam[:-4]+".n.sam"]
multiple_p=multiprocessing.Pool(2)
for file in p_n:
    ress=multiple_p.apply_async(run_ccs,args=(file,))
     
print('Waiting for all subprocesses done...')
multiple_p.close()
multiple_p.join()
print('All subprocesses done.')

