import os
from sys import argv
from collections import defaultdict, Counter


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




output1=open(out_path+"%s_new_major_geno_reads_filter_%s_v3.txt"%(file_prefix,umi_reads_filter),"w")
output1.write("Sample"+"\t"+"Ku_type"+"\t"+"Day"+"\t"+"Template"+"\t"+"Biorepeat"+"\t"
              +"umi3"+"\t"
              +"major_umi5"+"\t"
              +"major_geno"+"\t"
              +"syn_type"+"\t"
              +"reads"+"\n")
              
output2=open(out_path+"%s_new_m1_major_geno_reads_filter_%s_v3.txt"%(file_prefix,umi_reads_filter),"w")
output2.write("Sample"+"\t"+"Ku_type"+"\t"+"Day"+"\t"+"Template"+"\t"+"Biorepeat"+"\t"
              +"Genotype"+"\t"
              +"Syn_Type"+"\t"
              +"Major_umi_types_num"+"\t"
              +"Major_umi_types_num_sum"+"\n")

output3=open(out_path+"%s_new_m2_major_geno_reads_filter_%s_v3.txt"%(file_prefix,umi_reads_filter),"w")
output3.write("Sample"+"\t"+"Ku_type"+"\t"+"Day"+"\t"+"Template"+"\t"+"Biorepeat"+"\t"
              +"Genotype"+"\t"
              +"Syn_Type"+"\t"
              +"Major_umi_reads_num"+"\t"
              +"Major_umi_reads_num_sum"+"\n")


for filename in path_list:
    print (filename)
    
    with open(filename,"r") as ccs_file:
        sequences = []
        for line in ccs_file:
                line=line.strip()
                umi_g_list=line.split(" ")
                geno=umi_g_list[2]
                umi5=umi_g_list[3]#5'umi
                umi3=umi_g_list[4]#3'umi
                sequences.append((umi3, umi5, geno))


        reads_count = defaultdict(int)


        for umi3, umi5, geno in sequences:
            reads_count[(umi3, umi5, geno)] += 1


        data = [(umi3, umi5, geno, reads) for (umi3, umi5, geno), reads in reads_count.items()]


        umi3_reads = defaultdict(int)
        for umi3, umi5, geno, reads in data:
            umi3_reads[umi3] += reads


        filtered_umi3 = {umi3 for umi3, reads in umi3_reads.items() if reads > int(umi_reads_filter)}


        umi3_to_umi5 = defaultdict(Counter)
        for umi3, umi5, geno, reads in data:
            if umi3 in filtered_umi3:
                umi3_to_umi5[umi3][umi5] += reads

        umi3_major_umi5 = {}
        for umi3, umi5_counts in umi3_to_umi5.items():
            sorted_umi5 = sorted(umi5_counts.items(), key=lambda x: x[1], reverse=True)
            if len(sorted_umi5) > 1 and sorted_umi5[0][1] == sorted_umi5[1][1]:
                continue  
            umi3_major_umi5[umi3] = sorted_umi5[0][0]
            

        umi3_umi5_set = set((umi3, umi5) for umi3, umi5 in umi3_major_umi5.items())


        umi3_to_umi5_geno = defaultdict(lambda: defaultdict(Counter))


        for umi3, umi5, geno, reads in data:

            if (umi3, umi5) in umi3_umi5_set:
                umi3_to_umi5_geno[umi3][umi5][geno] += reads


        umi3_major_results = {}
        for umi3, umi5_data in umi3_to_umi5_geno.items():
            for umi5, geno_counts in umi5_data.items():
                
                if len(geno_counts) == 1:

                    geno = list(geno_counts.keys())[0]
                    umi3_major_results[umi3] = (umi5,geno)




        umi3_reads_gte_3_count = len(filtered_umi3)

        umi3_with_major_umi5_count = len(umi3_major_umi5)

        major_umi5_with_major_geno_count = len(umi3_major_results)

        umi5_counts = Counter(umi5 for umi3, (umi5, geno) in umi3_major_results.items())

        shared_umi5 = {umi5 for umi5, count in umi5_counts.items() if count > 1}

        umi3_with_shared_umi5 = [umi3 for umi3, (umi5, geno) in umi3_major_results.items() if umi5 in shared_umi5]



        geno_to_umi3 = defaultdict(set)
        for umi3, (umi5, geno) in umi3_major_results.items():
            geno_to_umi3[geno].add(umi3)

        geno_umi3_counts = {geno: len(umi3_set) for geno, umi3_set in geno_to_umi3.items()}


        umi3_major_results_set = set((umi3, umi5, geno) for umi3, (umi5, geno) in umi3_major_results.items())


        filtered_data = [
            (umi3, umi5, geno, reads)
            for umi3, umi5, geno, reads in data
            if (umi3, umi5, geno) in umi3_major_results_set
        ]


        geno_reads = defaultdict(int)


        for umi3, umi5, geno, reads in filtered_data:
            geno_reads[geno] += reads


        total_reads = sum(geno_reads.values())


        geno_sim_dict={}
        for g in geno_to_umi3.keys():
            genotype=""
            for i in range(len(wild_type)):
                if g[i]!= wild_type[i]:
                    s=wild_type[i]+str(i+1)+g[i]
                    genotype=genotype+s+" "
            if genotype=="":
                genotype="WT "
            geno_sim_dict[g]=genotype.strip()


        geno_mut_dict={}
        for geno_i in geno_to_umi3.keys():
            translate_mut=translate_dna(geno_i)
            if translate_mut!=translate_WT:
                geno_mut_dict[geno_i]="non_synonymous"
            else:
                geno_mut_dict[geno_i]="synonymous"


        for r in filtered_data:
            output1.write(filename.split(".")[0]+"\t"
                            +filename.split("_")[0]+"\t"
                            +filename.split("_")[1][1]+"\t"
                            +filename.split("_")[2].split("n")[0]+"\t"
                            +filename.split("_")[-1][0]+"\t"
                            +r[0]+"\t"
                            +r[1]+"\t"
                            +geno_sim_dict[(r[2])]+"\t"
                            +geno_mut_dict[(r[2])]+"\t"
                            +str(r[3])+"\n")


        for a in geno_umi3_counts.keys():
            output2.write(filename.split(".")[0]+"\t"
                            +filename.split("_")[0]+"\t"
                            +filename.split("_")[1][1]+"\t"
                            +filename.split("_")[2].split("n")[0]+"\t"
                            +filename.split("_")[-1][0]+"\t"
                            +geno_sim_dict[a]+"\t"
                            +geno_mut_dict[a]+"\t"
                            +str(geno_umi3_counts[a])+"\t"
                            +str(sum(geno_umi3_counts.values()))+"\n")


        for b in geno_reads.keys():
            output3.write(filename.split(".")[0]+"\t"
                            +filename.split("_")[0]+"\t"
                            +filename.split("_")[1][1]+"\t"
                            +filename.split("_")[2].split("n")[0]+"\t"
                            +filename.split("_")[-1][0]+"\t"
                            +geno_sim_dict[b]+"\t"
                            +geno_mut_dict[b]+"\t"
                            +str(geno_reads[b])+"\t"
                            +str(total_reads)+"\n")
            


output1.close()
output2.close()
output3.close()






