from sys import argv
import regex as re

raw_bam=argv[1]
print("raw_bam:",raw_bam)

mud_dict={"A":"T","T":"A","C":"G","G":"C"}
Reverse_complemrnt= lambda x:"".join([mud_dict[i] for i in x][::-1])

umi2_pattern=re.compile("ACTGGCCGTCGTTTTACA\w{20}GATCTCGGTG")
C_index_umi1_pattern=re.compile("CATACGAGAT\w{26}GGGTCCTTCTTCCTGAATACTC")
C_geno_pattern=re.compile("ATGCGGGTCATGGCGCCCCG[ATCG]+TATTGTGTTCTGTAGCCTGA")


outfile_C=open("%s_C_v4_pass3.final_result"%(raw_bam[:-4]),"w")

with open("%s.p_ccs_pass3.sam"%(raw_bam[:-4]),"r") as inputf:#positive strand ccs result      
    for posi in inputf:
        if posi[0]=="m":
            posi_line_list=posi.strip().split()
            posi_seq=posi_line_list[9]
            zmw_posi=posi_line_list[0].split("/")[1]
            umi_p2=re.search(umi2_pattern,posi_seq,flags=0)
            C_index_umi1=re.search(C_index_umi1_pattern,posi_seq,flags=0)
            C_geno=re.search(C_geno_pattern,posi_seq,flags=0)

            if C_geno is not None and C_index_umi1 is not None:
                if umi_p2 is not None:
                    C_index_p=C_index_umi1.group()[10:16]
                    C_umi_p11= C_index_umi1.group()[16:36]
                    C_geno1 = C_geno.group()[:]
                    C_umi_p22= umi_p2.group()[18:38]
                    outfile_C.write(">"+str(zmw_posi)+"p"+" "+C_index_p+" "+C_geno1+" "+C_umi_p11+" "+C_umi_p22+"\n")


                    
with open("%s.n_ccs_pass3.sam"%(raw_bam[:-4]),"r") as inputr:
    for nega in inputr:    
        if nega[0]=="m":
            nega_line_list=nega.strip().split()
            nega_seq=nega_line_list[9]
            zmw_nega=nega_line_list[0].split("/")[1]
            nega_reverse=Reverse_complemrnt(nega_seq)
            umi_n2=re.search(umi2_pattern,nega_reverse,flags=0)
            C_index2_umi1=re.search(C_index_umi1_pattern,nega_reverse,flags=0)
            C_geno2=re.search(C_geno_pattern,nega_reverse,flags=0)

            if C_geno2 is not None and C_index2_umi1 is not None:
                if  umi_n2 is not None:
                    C_index_n=C_index2_umi1.group()[10:16]
                    C_umi_n11 = C_index2_umi1.group()[16:36]
                    C_geno3 = C_geno2.group()[:]
                    C_umi_n22 = umi_n2.group()[18:38]
                    outfile_C.write(">"+str(zmw_nega)+"n"+" "+C_index_n+" "+C_geno3+" "+C_umi_n11+" "+C_umi_n22+"\n")

          
outfile_C.close()