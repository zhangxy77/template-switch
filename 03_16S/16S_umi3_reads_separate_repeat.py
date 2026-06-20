import os
from sys import argv

data_path,ku_type=argv[1:]

print("data_path:",data_path)
print("ku_type:",ku_type)

os.chdir(data_path)
os.chdir("..")
data_path_upDir=os.getcwd()
out_path=data_path_upDir+"/%s_Umi3_reads_sep/"%(ku_type)
os.mkdir(out_path)

os.chdir(data_path)
path_list=[f for f in os.listdir(data_path) if f.endswith('.txt')]
path_list.sort()

#output
output1=open(out_path+"%s_umi3_reads.txt"%(ku_type),"w")
output1.write("Sample"+"\t"+"Cycles"+"\t"+"Biorepeat"+"\t"+"Umi3"+"\t"+"reads"+"\n")

for filename in path_list:
    print (filename)

    with open(filename,"r") as ccs_file:
        umi1_umi2={}
        for line in ccs_file:
            line=line.strip()
            umi_g_list=line.split(" ")
            umi1=umi_g_list[4][:20]
            umi2=umi_g_list[3][:20]
            if umi1 in umi1_umi2.keys():
                umi1_umi2[umi1].append(umi2)
            else:
                umi1_umi2[umi1]=[umi2] 


        umi1_reads={}
        for u1 in umi1_umi2.keys():
            umi1_reads[u1]=len(umi1_umi2[u1])

        for u2 in umi1_umi2.keys():
            output1.write(filename[:-4]+"\t"
                          +filename.split("_")[1][:2]+"\t"
                          +filename.split("_")[0][-1]+"\t"
                          +str(u2)+"\t"
                          +str(umi1_reads[u2])+"\n"
                          )

output1.close()









