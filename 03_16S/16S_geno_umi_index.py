from sys import argv
import regex as re

raw_bam=argv[1]

mud_dict={"A":"T","T":"A","C":"G","G":"C"}
Reverse_complemrnt= lambda x:"".join([mud_dict[i] for i in x][::-1])

index_umi1_pattern=re.compile("CATACGAGAT\w{36}AG[AG]GTT[CT]GAT[CT][AC]TGGCTCAG")
geno_pattern=re.compile("AG[AG]GTT[CT]GAT[CT][AC]TGGCTCAG[ATCG]+AAGTCGTAACAAGGTA[GA]C[TC]")
umi2_pattern=re.compile("AAGTCGTAACAAGGTA[GA]C[TC]\w{20}GATCTCGGTG")

outfile=open("%s_v2.final_result"%(raw_bam[:-4]),"w")

with open("%s.p_ccs.sam"%(raw_bam[:-4]),"r") as inputf:
    for posi in inputf:
        if posi[0]=="m":
            posi_line_list=posi.strip().split()
            posi_seq=posi_line_list[9]
            zmw_posi=posi_line_list[0].split("/")[1]
            index_umi1=re.search(index_umi1_pattern,posi_seq,flags=0)
            geno=re.search(geno_pattern,posi_seq,flags=0)
            umi_p2=re.search(umi2_pattern,posi_seq,flags=0)
            if geno is not None and index_umi1 is not None:
                if umi_p2 is not None:
                    index_p=index_umi1.group()[10:26]
                    geno1 = geno.group()[:]
                    umi_p11= index_umi1.group()[26:46]
                    umi_p22= umi_p2.group()[19:39]
                    outfile.write(">"+str(zmw_posi)+"p"+" "+index_p+" "+geno1+" "+umi_p11+" "+umi_p22+"\n")
                    


with open("%s.n_ccs.sam"%(raw_bam[:-4]),"r") as inputr:
    for nega in inputr:    
        if nega[0]=="m":
            nega_line_list=nega.strip().split()
            nega_seq=nega_line_list[9]
            zmw_nega=nega_line_list[0].split("/")[1]
            nega_reverse=Reverse_complemrnt(nega_seq)
            index2_umi1=re.search(index_umi1_pattern,nega_reverse,flags=0)
            geno2=re.search(geno_pattern,nega_reverse,flags=0)
            umi_n2=re.search(umi2_pattern,nega_reverse,flags=0)
            if geno2 is not None and index2_umi1 is not None:
                if  umi_n2 is not None:
                    index_n=index2_umi1.group()[10:26]
                    geno3 = geno2.group()[:]
                    umi_n11 = index2_umi1.group()[26:46]
                    umi_n22 = umi_n2.group()[19:39]
                    outfile.write(">"+str(zmw_nega)+"n"+" "+index_n+" "+geno3+" "+umi_n11+" "+umi_n22+"\n")
                    
outfile.close()