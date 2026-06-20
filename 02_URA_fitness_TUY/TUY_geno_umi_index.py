from sys import argv
import regex as re

raw_bam=argv[1]

mud_dict={"A":"T","T":"A","C":"G","G":"C"}
Reverse_complemrnt= lambda x:"".join([mud_dict[i] for i in x][::-1])


umi1_pattern=re.compile("CATACGAGAT\w{20}TATTGATTGTAATTCTGTAA")

geno_pattern=re.compile("AACACACATAAACAAACAAA\w{804}ATAAGCGAATTTCTTATGAT")

umi2_index_pattern=re.compile("GAATTTCTTATGATTTATGA\w{26}GATCTCGGTG")

outfile=open("%s.final_result"%(raw_bam[:-4]),"w")
outfile2=open("%s.fasta"%(raw_bam[:-4]),"w")

with open("%s.p_ccs.sam"%(raw_bam[:-4]),"r") as inputf:
    for posi in inputf:
        if posi[0]=="m":
            posi_line_list=posi.strip().split()
            posi_seq=posi_line_list[9]
            zmw_posi=posi_line_list[0].split("/")[1]
            umi1_p=re.search(umi1_pattern,posi_seq,flags=0)
            geno_p=re.search(geno_pattern,posi_seq,flags=0)
            umi2_index_p=re.search(umi2_index_pattern,posi_seq,flags=0)
            if geno_p is not None and umi1_p is not None:
                if umi2_index_p is not None:
                    index_p=umi2_index_p.group()[40:46]
                    geno1 = geno_p.group()[20:824]
                    umi_p11= umi1_p.group()[10:30]
                    umi_p22= umi2_index_p.group()[20:40]
                    outfile.write(">"+str(zmw_posi)+"p"+" "+index_p+" "+geno1+" "+umi_p11+" "+umi_p22+"\n")
                    outfile2.write(">"+str(zmw_posi)+"p"+"\n"+geno1+"\n")


with open("%s.n_ccs.sam"%(raw_bam[:-4]),"r") as inputr:
    for nega in inputr:    
        if nega[0]=="m":
            nega_line_list=nega.strip().split()
            nega_seq=nega_line_list[9]
            zmw_nega=nega_line_list[0].split("/")[1]
            nega_reverse=Reverse_complemrnt(nega_seq)
            umi1_n=re.search(umi1_pattern,nega_reverse,flags=0)
            geno_n=re.search(geno_pattern,nega_reverse,flags=0)
            umi2_index_n=re.search(umi2_index_pattern,nega_reverse,flags=0)
            if geno_n is not None and umi1_n is not None:
                if  umi2_index_n is not None:
                    index_n=umi2_index_n.group()[40:46]
                    geno3 = geno_n.group()[20:824]
                    umi_n11 = umi1_n.group()[10:30]
                    umi_n22 = umi2_index_n.group()[20:40]
                    outfile.write(">"+str(zmw_nega)+"n"+" "+index_n+" "+geno3+" "+umi_n11+" "+umi_n22+"\n")
                    outfile2.write(">"+str(zmw_nega)+"n"+"\n"+geno3+"\n")
                   
                    

outfile.close()
outfile2.close()




