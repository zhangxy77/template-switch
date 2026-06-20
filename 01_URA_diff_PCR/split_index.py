import os
from sys import argv

result_file, index_info=argv[1:]

mud_dict={"A":"T","T":"A","C":"G","G":"C"}
Reverse_complemrnt= lambda x:"".join([mud_dict[i] for i in x][::-1])


out_path=os.getcwd()+"/splited_txt_file/"
os.mkdir(out_path)

out1=open("no_fitted_index_seq.txt","w")


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
        seq=seq_line.split()[1]
        seq_umi=seq_line.split()[2]
        seq_index=seq_line.split()[3]
        seq_index_rever=Reverse_complemrnt(seq_index)
        if seq_index_rever in index_dict.keys():
            filename=index_dict[seq_index_rever]
            with open(out_path+"%s.txt"%(filename),"a") as f1:
                f1.write(seq_name+"\n"+seq+" "+seq_umi+" "+seq_index+"\n")
                f1.close()
        else:
            out1.write(seq_line+"\n")


out1.close()



