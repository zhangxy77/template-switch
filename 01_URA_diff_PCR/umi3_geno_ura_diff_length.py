import os
from sys import argv
import numpy as np

sep_data_path,combine_data_path,outpath,ku_type,umi_reads_filter=argv[1:]

print("sep_data_path:",sep_data_path)
print("combine_data_path:",combine_data_path)
print("outpath:",outpath)
print("ku_type:",ku_type)
print("umi_reads_filter",umi_reads_filter)


out_path=outpath+"Umi3_geno_fig1_reads_%s/"%(umi_reads_filter)
os.mkdir(out_path)


output1=open(out_path+"%s_All_OtoM_ratio_3umi_%s.txt"%(ku_type,umi_reads_filter),"w")
output1.write("Sample"+"\t"+"Template"+"\t"+"Length"+"\t"+"Enzyme"+"\t"+"Cycles"+"\t"+"Biorepeat"+"\t"+"OTO1_umi_num"+"\t"+"OTOMore_umi_num"+"\t"+"Ratio"+"\n")

output2=open(out_path+"%s_minor_pro_3umi_%s.txt"%(ku_type,umi_reads_filter),"w")
output2.write("Sample"+"\t"+"Template"+"\t"+"Length"+"\t"+"Enzyme"+"\t"+"Cycles"+"\t"+"Biorepeat"+"\t"+"Umi3_reads"+"\t"+"Major_geno_reads"+"\t"+"Minor_geno_reads"+"\t"+"Minor_proportion"+"\n")

output3=open(out_path+"Com_repeat_%s_3umi_geno_types_3umi_%s.txt"%(ku_type,umi_reads_filter),"w")
output3.write("Sample"+"\t"+"Template"+"\t"+"Length"+"\t"+"Enzyme"+"\t"+"Cycles"+"\t"+"Umi3"+"\t"+"Geno_type_num"+"\n")

output4=open(out_path+"Com_repeat_%s_3umi_geno_major_pro_3umi_%s.txt"%(ku_type,umi_reads_filter),"w")
output4.write("Sample"+"\t"+"Template"+"\t"+"Length"+"\t"+"Enzyme"+"\t"+"Cycles"+"\t"+"Umi3"+"\t"+"Umi3_reads"+"\t"+"Major_geno_reads"+"\t"+"Minor_geno_reads"+"\t"+"Major_proportion"+"\n")



os.chdir(sep_data_path)
sep_path_list=[f for f in os.listdir(sep_data_path) if f.endswith('.txt')]
sep_path_list.sort()


for filename in sep_path_list:
    print (filename)

    with open(filename,"r") as ccs_file:
        umi1_umi2={}
        for line in ccs_file:
            line=line.strip()
            if line[0]=="n" or line[0]=="p" :
                umi_g_list=line.split(" ")
                umi1=umi_g_list[1][:20]
                umi2=umi_g_list[0][2::]
                if umi1 in umi1_umi2.keys():
                    umi1_umi2[umi1].append(umi2)
                else:
                    umi1_umi2[umi1]=[umi2]
            else:
                if line[0] !=">":
                    print(line)




        umi1_reads={}
        umi1_reads_filter_umi2={}
        for u1 in umi1_umi2.keys():
            if len(umi1_umi2[u1])>int(umi_reads_filter):
                umi1_reads[u1]=len(umi1_umi2[u1])
                umi1_reads_filter_umi2[u1]=umi1_umi2[u1]



        oTo=0
        oTm=0

        major_reads=0
        minor_reads=0
        major_umi_reads=0


        for u3 in umi1_reads_filter_umi2.keys():
            if len(np.unique(umi1_reads_filter_umi2[u3]))==1:
                oTo=oTo+1
                major_umi_reads=major_umi_reads+umi1_reads[u3]
                major_reads=major_reads+umi1_reads[u3]

            else:
                oTm=oTm+1

                umi2_count_dict={}
                uni_umi2=np.unique(umi1_umi2[u3])#
                for umi2 in uni_umi2:
                    umi2_count_dict[umi2]=umi1_umi2[u3].count(umi2)
                sort_umi2_dict=sorted(umi2_count_dict.items(),key=lambda umi2_count_dict:umi2_count_dict[1],reverse=True)
                major_umi2_reads=sort_umi2_dict[0][1]
                second_umi2_reads=sort_umi2_dict[1][1]
                rest_geno_reads=0
                if major_umi2_reads>second_umi2_reads:
                    major_reads=major_reads+major_umi2_reads
                    major_umi_reads=major_umi_reads+umi1_reads[u3]
                    for rest_geno in sort_umi2_dict[1:]:
                        rest_geno_reads+=rest_geno[1]
                    minor_reads=minor_reads+rest_geno_reads


        ratio=oTm/oTo
        minor_perc=minor_reads/major_umi_reads


        samplename = filename.replace("_", "-")



        output1.write(samplename[:-4]+"\t"
                     +samplename.split("-")[0]+"\t"
                     +samplename.split("-")[1]+"\t"
                     +samplename.split("-")[2]+"\t"
                     +samplename.split("-")[3]+"\t"
                     +samplename.split("-")[4][:-4]+"\t"
                      +str(oTo)+"\t"
                      +str(oTm)+"\t"
                      +str(ratio)+"\n")

        output2.write(samplename[:-4]+"\t"
                     +samplename.split("-")[0]+"\t"
                     +samplename.split("-")[1]+"\t"
                     +samplename.split("-")[2]+"\t"
                     +samplename.split("-")[3]+"\t"
                     +samplename.split("-")[4][:-4]+"\t"
                      +str(major_umi_reads)+"\t"
                      +str(major_reads)+"\t"
                      +str(minor_reads)+"\t"
                      +str(minor_perc)+"\n")

output1.close()
output2.close()

print("Sep data processing finished, now start combine data processing...")


os.chdir(combine_data_path)
combine_path_list=[f for f in os.listdir(combine_data_path) if f.endswith('.txt')]
combine_path_list.sort()

for filename in combine_path_list:
    print (filename)


    with open(filename,"r") as ccs_file:
        umi1_umi2={}
        for line in ccs_file:
            line=line.strip()
            if line[0]=="n" or line[0]=="p" :
                umi_g_list=line.split(" ")
                umi1=umi_g_list[1][:20]
                umi2=umi_g_list[0][2::]
                if umi1 in umi1_umi2.keys():
                    umi1_umi2[umi1].append(umi2)
                else:
                    umi1_umi2[umi1]=[umi2] 
            else:
                if line[0] !=">":
                    print(line)


        umi1_reads={}
        umi1_reads_filter_umi2={}
        umi1_types={}
        for u1 in umi1_umi2.keys():
            if len(umi1_umi2[u1])>int(umi_reads_filter):
                umi1_reads[u1]=len(umi1_umi2[u1])
                umi1_reads_filter_umi2[u1]=umi1_umi2[u1]
                umi1_types[u1]=len(np.unique(umi1_umi2[u1]))

        umi1_major_umi2={}
        umi1_major_umi2_reads={}
        umi1_minor_umi2_reads={}
        for u4 in umi1_types.keys():
            if umi1_types[u4]==1:
                umi1_major_umi2[u4]=1
                umi1_major_umi2_reads[u4]=umi1_reads[u4]
                umi1_minor_umi2_reads[u4]=0
            else:
                umi2_count_dict={}
                uni_umi2=np.unique(umi1_umi2[u4])
                for umi2 in uni_umi2:
                    umi2_count_dict[umi2]=umi1_umi2[u4].count(umi2)
                sort_umi2_dict=sorted(umi2_count_dict.items(),key=lambda umi2_count_dict:umi2_count_dict[1],reverse=True)
                major_umi2_reads=sort_umi2_dict[0][1]
                second_umi2_reads=sort_umi2_dict[1][1]
                rest_geno_reads=0
                if major_umi2_reads>second_umi2_reads:
                    pro=major_umi2_reads/umi1_reads[u4]
                    umi1_major_umi2[u4]=pro
                    umi1_major_umi2_reads[u4]=major_umi2_reads
                    for rest_geno in sort_umi2_dict[1:]:
                        rest_geno_reads+=rest_geno[1]
                    umi1_minor_umi2_reads[u4]=rest_geno_reads

        samplename = filename.replace("_", "-")


        for u5 in umi1_types.keys():
            output3.write(samplename[:-4]+"\t"
                          +samplename.split("-")[0]+"\t"
                          +samplename.split("-")[1]+"\t"
                          +samplename.split("-")[2]+"\t"
                          +samplename.split("-")[3][:-4]+"\t"
                          +str(u5)+"\t"
                          +str(umi1_types[u5])+"\n")

        for u6 in umi1_major_umi2.keys():
            output4.write(samplename[:-4]+"\t"
                          +samplename.split("-")[0]+"\t"
                          +samplename.split("-")[1]+"\t"
                          +samplename.split("-")[2]+"\t"
                          +samplename.split("-")[3][:-4]+"\t"
                          +str(u6)+"\t"
                          +str(umi1_reads[u6])+"\t"
                          +str(umi1_major_umi2_reads[u6])+"\t"
                          +str(umi1_minor_umi2_reads[u6])+"\t"
                          +str(umi1_major_umi2[u6])+"\n")


            

output3.close()
output4.close()         

print("Combine data processing finished, all done!")

