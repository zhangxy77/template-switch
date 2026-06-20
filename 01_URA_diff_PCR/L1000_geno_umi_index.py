from sys import argv
import regex as re

raw_bam=argv[1]

mud_dict={"A":"T","T":"A","C":"G","G":"C"}
Reverse_complemrnt= lambda x:"".join([mud_dict[i] for i in x][::-1])


geno_pattern=re.compile("TGCAAGTCCGGTTGC\w{751}ATAAGCGAATTTCTTATGAT")
umi_pattern=re.compile("TCGCTCTTATTGACCACACC\w{25}GCTTCGGCAGCACATATACT")
index_pattern=re.compile("GGCAGCACATATACTAAGAT\w{16}GAGATAC")


outfile=open("L1000_geno_umi_index_version2_%s.final_result"%(raw_bam[:-4]),"w")
outfile2=open("L1000_geno_umi_index_version2_%s.fasta"%(raw_bam[:-4]),"w")
with open("%s.n_ccs.sam"%(raw_bam[:-4]),"r") as inputr:
    for nega in inputr:    
        if nega[0]=="m":
            nega_line_list=nega.strip().split()
            nega_seq=nega_line_list[9]
            zmw_nega=nega_line_list[0].split("/")[1]
            nega_reverse=Reverse_complemrnt(nega_seq)
            geno2=re.search(geno_pattern,nega_reverse,flags=0)
            umi2=re.search(umi_pattern,nega_reverse,flags=0)
            index2=re.search(index_pattern,nega_reverse,flags=0)
            if geno2 is not None and umi2 is not None:
                if index2 is not None:
                    geno3 = geno2.group()[15:766]
                    umi3 = umi2.group()[20:40]
                    index3=index2.group()[20:36]
                    outfile.write(">"+str(zmw_nega)+"n"+" "+"n_"+geno3+" "+umi3+" "+index3+"\n")
                    outfile2.write(">"+str(zmw_nega)+"n"+"\n"+geno3+"\n")
with open("%s.p_ccs.sam"%(raw_bam[:-4]),"r") as inputf:
    for posi in inputf:
        if posi[0]=="m":
            posi_line_list=posi.strip().split()
            posi_seq=posi_line_list[9]
            zmw_posi=posi_line_list[0].split("/")[1]
            geno=re.search(geno_pattern,posi_seq,flags=0)
            umi=re.search(umi_pattern,posi_seq,flags=0)
            index=re.search(index_pattern,posi_seq,flags=0)
            if geno is not None and umi is not None:
                if index is not None:#index
                    geno1 = geno.group()[15:766]
                    umi1 = umi.group()[20:40]
                    index1 = index.group()[20:36]
                    outfile.write(">"+str(zmw_posi)+"p"+" "+"p_"+geno1+" "+umi1+" "+index1+"\n")
                    outfile2.write(">"+str(zmw_posi)+"p"+"\n"+geno1+"\n")
                

outfile.close()
outfile2.close()