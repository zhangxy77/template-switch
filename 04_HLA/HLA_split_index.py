import os
from sys import argv

result_file, index_info,geno_name,out_path=argv[1:]

mud_dict={"A":"T","T":"A","C":"G","G":"C"}
Reverse_complemrnt= lambda x:"".join([mud_dict[i] for i in x][::-1])

print("result_file:",result_file)
print("index_info:",index_info)
print("geno_name:",geno_name)
print("out_path:",out_path)

out1=open("%s_no_fitted_index_seq.txt"%(geno_name),"w")

index_dict={}
with open(index_info,"r") as index_file:
    for index in index_file:
        index=index.strip()
        index_name=index.split()[0]
        index_seq=index.split()[1]
        index_dict[index_seq]=index_name

with open(result_file,"r") as seq_file:
    for seq_line in seq_file:
        seq_line=seq_line.strip()
        seq_name=seq_line.split()[0]
        seq_index=seq_line.split()[1]
        seq=seq_line.split()[2]
        seq_umi1=seq_line.split()[3]
        seq_umi2=seq_line.split()[4]
        if seq_index in index_dict.keys():
            filename=index_dict[seq_index]
            with open(out_path+"%s.txt"%(filename),"a") as f1:
                f1.write(seq_name+" "+seq_index+" "+seq+" "+seq_umi1+" "+seq_umi2+"\n")
                f1.close()
            with open(out_path+"%s.fasta"%(filename),"a") as f2:
                f2.write(seq_name+"\n"+seq+"\n")
                f2.close()
        else:
            out1.write(seq_line+"\n")

out1.close()