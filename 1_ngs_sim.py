import os
import random
import string

#function to generate simulation reads
#parameters are genome fasta file, length of reads, number of reads to be generated
def gen_sim_reads(ref_gen_fa, read_len, n_reads):
    ref_seq = ''                                    #ref_seq as string
    with open(ref_gen_fa) as f:
        for lines in f:
            lines = str(lines).strip('\n').upper()  #preparing the sequence
            if('>' not in lines):                   #ignoring id line
                ref_seq = ref_seq+str(lines)        #ref_seq string

    ref_seq = ref_seq.replace('N', '')              #removing N nucleotide

    ref_seq_length = len(ref_seq)
    #list to store generated reads
    sim_reads = []  
    for i in range(n_reads):
        rand_start_ind = random.randint(0, ref_seq_length-read_len-1)   #generating random int
        sim_reads.append(ref_seq[rand_start_ind : rand_start_ind+read_len])
    return(sim_reads)

def add_rand(list_reads, perc):
    n_reads = len(list_reads)                         # size of the list generated(100000)
    read_len = len(list_reads[0])                     # size of each string(element) in the list(50)
    total_reads_length = n_reads*read_len             # total number of characters
    total_rand = total_reads_length*perc/100          # total number of errors to be made

    pos_visited = []                                  # maintains list indexes already visited
    count = 0
    bases = ['A', 'T', 'G', 'C']
    while(count < total_rand):  
        i = random.randint(0,total_reads_length)      # generating random number between(0-5000000)
        if i in pos_visited: 
            pass
        else:
            pos_visited.append(i)
            bases_temp = bases[:]

            read_num = (int)(i/read_len)               # determine the index in the list where changes are to be made
            
            read_pos=50
            read_pos = i%read_pos                       # determine the position of the character in the list[read_num] where changes are to be made            # print (read_num,read_pos)
            base_flip = list_reads[read_num][read_pos] 
            bases_temp.remove(base_flip)                # remove the char from the bases_temp that was present in the loc where change is to be made
            base_replace = random.choice(bases_temp)    # chose the new char 
            tmp=list(list_reads[read_num])              # convert the string at that index in list to a list 
            tmp[read_pos]=base_replace                  # replace the old char by a new random char 
            list_reads[read_num]="".join(tmp)           # modify the string 
            count+=1 
    return(list_reads)  


def random_char(y):
    return ''.join(random.choice(string.ascii_letters) for x in range(y))

print ("***************Reading the fasta file***************")
print ("***************Takes a while to read hg38.fa file ***************")

reads = gen_sim_reads('hg38_truncated.fa', 50, 100000) 

print ("***************Successfully read the file***************")
list_reads = "\n".join(reads)
# print (list_reads)
reads = add_rand(reads, 1)                              # passing reads and error percentage
print ("********************Performed random replacements****************")
for i in range(len(reads)):
     print (i,reads[i])

last_line={0:"!",1:"\"",2:"#",3:"$",4:"%",5:"&",6:"\'",7:"(",8:")",9:"*",10:"+",11:",",12:"-",13:".",14:"/",15:"0",16:"1",17:"2",18:"3",19:"4",20:"5",21:"6",22:"7",23:"8",24:"9",25:":",26:";",27:"<",28:"=",29:"?",30:"?",31:"@",32:"A",33:"B",34:"C",35:"D",36:"E",37:"F",38:"G",39:"H",40:"I",41:"J"}
YorN={0:"Y",1:"N"}                                      #YorN is created for the 3rd last character of the 1st line (4 lines wala)

final_corrected_write=open("fastqfile.txt","w")

for i in range(len(reads)):
    # converting fasta file to fastq
    line1 = "@" + random_char(3).upper() + ":" + str(random.randint(0,9)) + ":" + random_char(3).upper() + str(random.randint(0,9)) + ":" 
    line1 = line1 + str(random.randint(10,99)) + ":" + str(random.randint(1000,9999)) + " " + str(random.randint(1000,9999)) +":" + str(random.randint(0,9)) + ":" 
    line1 = line1 + str(YorN[random.randint(0,1)]) + ":" + str(random.randint(0,9)) + ":" + str(random.randint(0,9)) + "\n"
    
    final_corrected_write.write(line1)
    line2=reads[i] + "\n"
    final_corrected_write.write(line2)
    final_corrected_write.write("+\n")
    line4=""
    for j in range(0,50):
        line4=line4+last_line[random.randint(0,41)]
    line4=line4+"\n"
    final_corrected_write.write(line4)

print ("***************Successfully created the final file***************")
print ("***************fastqfile.txt has the fastq file form style***************")

final_corrected_write.close()
