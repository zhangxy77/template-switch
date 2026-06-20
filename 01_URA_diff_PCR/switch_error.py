import os
from sys import argv
import numpy as np

data_path,leng,umi_reads_filter=argv[1:]
print("data_path:",data_path)
print("length:",leng)
print("umi_reads_filter:",int(umi_reads_filter))

wild_type="ATCGAAACTTTTACGTTTAATGGATGAAAAGAAGACCAATTTGTGTGCTTCTCTTGACGTTCGTTCGACTGATGAGCTATTGAAACTTGTTGAAACGTTGGGTCCATACATTTGCCTTTTGAAAACACACGTTGATATCTTGGATGATTTCAGTTATGAGGGTACTGTCGTTCCATTGAAAGCATTGGCAGAGAAATACAAGTTCTTGATATTTGAGGACAGAAAATTCGCCGATATCGGTAACACAGTCAAATTACAATATACATCGGGCGTTTACCGTATCGCAGAATGGTCTGATATCACCAACGCCCACGGGGTTACTGGTGCTGGTATTGTTGCTGGCTTGAAACAAGGTGCGCAAGAGGTCACCAAAGAACCAAGGGGATTATTGATGCTTGCTGAATTGTCTTCCAAGGGTTCTCTAGCACACGGTGAATATACTAAGGGTACCGTTGATATTGCAAAGAGTGATAAAGATTTCGTTATTGGGTTCATTGCTCAGAACGATATGGGAGGAAGAGAAGAAGGGTTTGATTGGCTAATCATGACCCCAGGTGTAGGTTTAGACGACAAAGGCGATGCATTGGGTCAGCAGTACAGAACCGTCGACGAAGTTGTAAGTGGTGGATCAGATATCATCATTGTTGGCAGAGGACTTTTCGCCAAGGGTAGAGATCCTAAGGTTGAAGGTGAAAGATACAGAAATGCTGGATGGGAAGCGTACCAAAAGAGAATCAGCGCTCCCCATTAA"
sp_ic="-"

print(wild_type)


def get_keys_from_value(dictionary, value):
    return [k for k, v in dictionary.items() if v == value]

os.chdir(data_path)
os.chdir("..")
data_path_upDir=os.getcwd()
out_path=data_path_upDir+"/combine_repeat_umi_switch_error_use_freq_v5_umi_2/"
os.mkdir(out_path)

os.chdir(data_path)
path_list=os.listdir(data_path)
path_list.sort()


output2=open(out_path+"%s_ComRepeat_switch_error_v5_umi_%s.txt"%(leng,umi_reads_filter),"w")
output2.write("Sample"+"\t"+"Template"+"\t"+"Length"+"\t"+"Enzyme"+"\t"+"Cycles"+"\t"+
              "umi"+"\t"+
              "major_loc"+"\t"+
              "switch_rate_and_error_rate"+"\t"+
              "error_rate"+"\t"+
              "switch_rate"+"\t"+
              "umi_reads"+"\t"+
              "major_reads"+"\t"+
              "before_1major_type"+"\t"+
              "before_1major_reads"+"\t"+
              "after_major1_type"+"\t"+
              "after_major1_reads"+"\t"+
              "before_1major1_after_type"+"\t"+
              "before_1major1_after_reads"+"\t"+
              "other_type"+"\t"+
              "other_type_reads"+"\n")



for filename in path_list:
    print (filename)

    with open(filename,"r") as ccs_file:
        umi_genot_dict={}
        for line in ccs_file:
            line=line.strip()
            if line[0]=="n" or line[0]=="p" :
                umi_g_list=line.split(" ")
                umi_pre=umi_g_list[1][:20]
                geno_pre=umi_g_list[0][2::]
                if umi_pre in umi_genot_dict.keys():
                    umi_genot_dict[umi_pre].append(geno_pre)
                else:
                    umi_genot_dict[umi_pre]=[geno_pre] 
            else:
                if line[0] !=">":
                    print(line)


        new_umi_genot_dict={}

        for umi in umi_genot_dict.keys():
            for geno_i in umi_genot_dict[umi][0:]:
                genotype=""
                for g in range(len(wild_type)):
                    if geno_i[g]!= wild_type[g]:
                        s=wild_type[g]+str(g+1)+geno_i[g]
                        genotype=genotype+s+" "
                if genotype=="":
                    genotype="WT "
                if umi in new_umi_genot_dict.keys():
                    new_umi_genot_dict[umi].append(genotype)
                else:
                    new_umi_genot_dict[umi]=[genotype]


        oneToone_umi_list=[]
        More_umi_list=[]
        umi_mut_num_dict={}
        umi_major_dict={}
        umi_nomajor_list=[]

        for umi2 in new_umi_genot_dict.keys():
            if len(np.unique(new_umi_genot_dict[umi2]))==1:
                only_geno=new_umi_genot_dict[umi2][0]
                if new_umi_genot_dict[umi2].count(only_geno)>int(umi_reads_filter):
                    oneToone_umi_list.append(umi2)
            else:
                if len(new_umi_genot_dict[umi2])>int(umi_reads_filter):
                    More_umi_list.append(umi2)
                    geno_count_dict={}
                    uni_gen=np.unique(new_umi_genot_dict[umi2])
                    for g2 in uni_gen:
                        geno_count_dict[g2]=new_umi_genot_dict[umi2].count(g2)
                    sort_geno_dict=sorted(geno_count_dict.items(),key=lambda geno_count_dict:geno_count_dict[1],reverse=True)
                    major_reads=sort_geno_dict[0][1]
                    second_reads=sort_geno_dict[1][1]
                    if major_reads>second_reads:
                        major=sort_geno_dict[0][0].strip()
                        umi_major_dict[umi2]=major
                        if major =="WT":
                            mut_num=0
                            umi_mut_num_dict[umi2]=mut_num
                        else:
                            mut_num=len(major.strip().split(" "))
                            umi_mut_num_dict[umi2]=mut_num
                    else:
                        umi_nomajor_list.append(umi2)                    


        umi_major_1=get_keys_from_value(umi_mut_num_dict, 1)
        umi_have_1major={}
        umi_no_1major={}
        umi_special_minor={}
        umi1_s_e={}

        for umi3 in umi_major_1[:]:
            umi_reads=len(new_umi_genot_dict[umi3])
            major1=umi_major_dict[umi3]+" "
            major1_reads=new_umi_genot_dict[umi3].count(major1)
            all_geno=np.unique(new_umi_genot_dict[umi3])
            minor_list = [x for x in all_geno if x != major1]


            umi_special_minor[umi3]=[]

            s=0 
            e=0 
            ss=0
            


            s1=0
            s1_reads=0
            e1=0
            e1_reads=0
            c1=0
            c1_reads=0
            c2=0
            c2_reads=0


            for m in minor_list[:]:
                m_list=m.strip().split(" ")
                if major1.strip() in m_list:
                    major_loc=int(major1.strip()[1:-1])
                    minor_mut_list=[x for x in m_list if x != major1.strip()]
                    minor_mut_list_len=len(minor_mut_list)

                    if minor_mut_list_len==1:
                        for minor_mut in minor_mut_list[:]:
                            loc=int(minor_mut[1:-1])
                            reads=new_umi_genot_dict[umi3].count(m)
                            if loc < major_loc :
                                s=s+((reads/umi_reads)/(major_loc-1))
                                s1=s1+1
                                s1_reads=s1_reads+reads
                            elif loc > major_loc :
                                e=e+((reads/umi_reads)/(len(wild_type)-major_loc))
                                e1=e1+1
                                e1_reads=e1_reads+reads
                            else:
                                umi_special_minor[umi3].append(m)

                    elif minor_mut_list_len==2:
                        loc1=int(minor_mut_list[0][1:-1])
                        loc2=int(minor_mut_list[1][1:-1])
                        reads=new_umi_genot_dict[umi3].count(m)
                        if loc1<major_loc and loc2>major_loc:
                            s=s+((reads/umi_reads)/(major_loc-1))
                            e=e+((reads/umi_reads)/(len(wild_type)-major_loc))
                            c1=c1+1
                            c1_reads=c1_reads+reads
                        elif loc2<major_loc and loc1>major_loc:
                            s=s+((reads/umi_reads)/(major_loc-1))
                            e=e+((reads/umi_reads)/(len(wild_type)-major_loc))
                            c1=c1+1
                            c1_reads=c1_reads+reads
                        else:
                            c2=c2+1
                            c2_reads=c2_reads+new_umi_genot_dict[umi3].count(m)         

                    else:
                        c2=c2+1
                        c2_reads=c2_reads+new_umi_genot_dict[umi3].count(m)
                    
                    
                    if umi3 in umi_have_1major.keys():
                        umi_have_1major[umi3].append(m)
                    else:
                        umi_have_1major[umi3]=[m]


                else:
                    if umi3 in umi_no_1major.keys():
                        umi_no_1major[umi3].append(m)
                    else:
                        umi_no_1major[umi3]=[m]
            
            ss=s-e
            
            
            if umi3 in umi_have_1major.keys():
                umi1_s_e[umi3]=[major_loc,s,e,ss,umi_reads,major1_reads,s1,s1_reads,e1,e1_reads,c1,c1_reads,c2,c2_reads]

    
        
        for umi4 in umi1_s_e.keys():
            output2.write(filename[:-4]+"\t"+
                          filename.split(sp_ic)[0]+"\t"+
                          filename.split(sp_ic)[1]+"\t"+
                          filename.split(sp_ic)[2]+"\t"+
                          filename.split(sp_ic)[3][:-4]+"\t"+
                          str(umi4)+"\t"+
                          str(umi1_s_e[umi4][0])+"\t"+
                          str(umi1_s_e[umi4][1])+"\t"+
                          str(umi1_s_e[umi4][2])+"\t"+
                          str(umi1_s_e[umi4][3])+"\t"+
                          str(umi1_s_e[umi4][4])+"\t"+
                          str(umi1_s_e[umi4][5])+"\t"+
                          str(umi1_s_e[umi4][6])+"\t"+
                          str(umi1_s_e[umi4][7])+"\t"+
                          str(umi1_s_e[umi4][8])+"\t"+
                          str(umi1_s_e[umi4][9])+"\t"+
                          str(umi1_s_e[umi4][10])+"\t"+
                          str(umi1_s_e[umi4][11])+"\t"+
                          str(umi1_s_e[umi4][12])+"\t"+
                          str(umi1_s_e[umi4][13])+"\n")


output2.close()



















