import os
from sys import argv
from collections import Counter

data_path,umi_reads_filter,out_path,file_prefix=argv[1:]
print("data_path:",data_path)
print("umi_reads_filter:",int(umi_reads_filter))
print("out_path:",out_path)
print("file_prefix:",file_prefix)

wild_type="ATGTCCACAAAATCATATACCAGTAGAGCTGAGACTCATGCAAGTCCGGTTGCATCGAAACTTTTACGTTTAATGGATGAAAAGAAGACCAATTTGTGTGCTTCTCTTGACGTTCGTTCGACTGATGAGCTATTGAAACTTGTTGAAACGTTGGGTCCATACATTTGCCTTTTGAAAACACACGTTGATATCTTGGATGATTTCAGTTATGAGGGTACTGTCGTTCCATTGAAAGCATTGGCAGAGAAATACAAGTTCTTGATATTTGAGGACAGAAAATTCGCCGATATCGGTAACACAGTCAAATTACAATATACATCGGGCGTTTACCGTATCGCAGAATGGTCTGATATCACCAACGCCCACGGGGTTACTGGTGCTGGTATTGTTGCTGGCTTGAAACAAGGTGCGCAAGAGGTCACCAAAGAACCAAGGGGATTATTGATGCTTGCTGAATTGTCTTCCAAGGGTTCTCTAGCACACGGTGAATATACTAAGGGTACCGTTGATATTGCAAAGAGTGATAAAGATTTCGTTATTGGGTTCATTGCTCAGAACGATATGGGAGGAAGAGAAGAAGGGTTTGATTGGCTAATCATGACCCCAGGTGTAGGTTTAGACGACAAAGGCGATGCATTGGGTCAGCAGTACAGAACCGTCGACGAAGTTGTAAGTGGTGGATCAGATATCATCATTGTTGGCAGAGGACTTTTCGCCAAGGGTAGAGATCCTAAGGTTGAAGGTGAAAGATACAGAAATGCTGGATGGGAAGCGTACCAAAAGAGAATCAGCGCTCCCCATTAA"


gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}


def translate_dna(dna):
    amino_acid_sequence = ""
    for start in range(0,len(dna) - 2, 3):
        stop = start + 3
        codon = dna[start:stop]
        aa = gencode.get(codon.upper(),'X') 
        amino_acid_sequence = amino_acid_sequence + aa
    return(amino_acid_sequence)

translate_WT=translate_dna(wild_type)

os.chdir(data_path)
path_list=os.listdir(data_path)
path_list.sort()

output1=open(out_path+"%s_Fitness_m2_pair_reads_filter_%s.txt"%(file_prefix,umi_reads_filter),"w")
output1.write("Sample"+"\t"+"Ku_type"+"\t"+"Day"+"\t"+"Template"+"\t"+"Biorepeat"+"\t"
              +"Genotype"+"\t"
              +"Syn_Type"+"\t"
              +"Major_umi_reads_num"+"\t"
              +"Major_umi_reads_num_sum"+"\n")


for filename in path_list:
    print(filename)

    with open(filename,"r") as ccs_file:
        geno_list=[]
        for line in ccs_file:
            line=line.strip()
            line_list=line.split(" ")
            geno_pre=line_list[2]
            geno_list.append(geno_pre)
        

        geno_reads=Counter(geno_list)
        new_geno_list=[geno for geno in geno_list if geno_reads[geno]>int(umi_reads_filter)]
        

        new_geno_reads=Counter(new_geno_list)

        reads_sum=len(new_geno_list)


        geno_sim_dict={}
        for g in new_geno_reads.keys():
            genotype=""
            for i in range(len(wild_type)):
                if g[i]!= wild_type[i]:
                    s=wild_type[i]+str(i+1)+g[i]
                    genotype=genotype+s+" "
            if genotype=="":
                genotype="WT "
            geno_sim_dict[g]=genotype.strip()


        geno_mut_dict={}
        for geno_i in new_geno_reads.keys():
            translate_mut=translate_dna(geno_i)
            if translate_mut!=translate_WT:
                geno_mut_dict[geno_i]="non_synonymous"
            else:
                geno_mut_dict[geno_i]="synonymous"


        for geno in new_geno_reads.keys():
            output1.write(filename.split(".")[0]+"\t"
                          +filename.split("_")[0]+"\t"
                          +filename.split("_")[1][1]+"\t"
                          +filename.split("_")[2].split("n")[0]+"\t"
                          +filename.split("_")[-1][0]+"\t"
                          +geno_sim_dict[geno]+"\t"
                          +geno_mut_dict[geno]+"\t"
                          +str(new_geno_reads[geno])+"\t"
                          +str(reads_sum)+"\n")


output1.close()






