import os
from sys import argv
from itertools import groupby

needle_file,matchtype,seq_file=argv[1:]


print("needle_file:",needle_file)
print("matchtype:",matchtype)
print("seq_file:",seq_file)

output1=open("no_indel.txt","w")
output2=open("indel.txt","w")

def split_text(s):
    for k, g in groupby(s,str.isalpha):

        yield ''.join(g)

with open(needle_file,"r") as needle_file1:
    for line in needle_file1:
        line=line.strip()
        check=0
        match_info=line.split()[2]
        match_info_split=list(split_text(match_info))
        for list_index, list_element in enumerate(match_info_split):
            if list_element == matchtype[:-1]:
                list_index1 = list_index
                aa=match_info_split[list_index1]+match_info_split[list_index1+1]
                if aa == matchtype:
                    check=1
                    output1.write(line+"\n")
        if check == 0:
            output2.write(line+"\n")

output1.close()
output2.close()

os.system('''time awk 'NR==FNR{a[">"$1]; next} $1 in a {print $0}' %s %s >> no_indel_seq.txt'''%(str(output1.name),seq_file))

