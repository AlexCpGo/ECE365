import numpy as np
from collections import OrderedDict

class Lab1(object):
    def parse_reads_illumina(self,reads) :
        '''
        Input - Illumina reads file as a string
        Output - list of DNA reads
        '''
        #start code here
        dna_reads_illumina=[]
        l=0 # l stands for the line number of the reads
        begin=0
        end=0
        beginlist=[]
        endlist=[]
        for i in range(len(reads)):
            if reads[i] == '\n':
                l+=1
                if l % 4 == 1:
                    begin=i+1
              
                elif l % 4 == 2:
                    end=i              
                if end>begin and reads[begin:end] not in dna_reads_illumina:
                    beginlist.append(begin)
                    endlist.append(end)
        for i in range(len(beginlist)):
            if i % 3 == 0:
                dna_reads_illumina.append(reads[beginlist[i]:endlist[i]])
        return dna_reads_illumina
        #end code here

    def unique_lengths(self,dna_reads) :
        '''
        Input - list of dna reads
        Output - set of counts of reads
        '''
        #start code here
        DNAreads_length = set()
        for i in dna_reads:
            if len(i) not in DNAreads_length:
                DNAreads_length.add(len(i))
        return DNAreads_length
        #end code here

    def check_impurity(self,dna_reads) :
        '''
        Input - list of dna reads
        Output - list of reads which have impurities, a set of impure chars 
        '''
        #start code here
        flag=0
        # Check which form is the dna_reads(capital/lower)
        for char in dna_reads:
            if char.upper()!=char:
                change=1
           
        Pure_list={'A','C','T','G'}
        if flag ==1:
            Pure_list={'a','c','t','g'}
        Impure_list=[]
        Impure_chars={'\n'}
        for item in dna_reads:
            set_items={item[i] for i in range(len(item))}
            set_extra=set_items - Pure_list
            length=len(set_extra)
            
            if length!=0:
                Impure_list.append(item)  
            for i in set_extra:   
                Impure_chars.add(i)
               
        Impure_chars.remove('\n')
        return Impure_list,Impure_chars
        #end code here

    def get_read_counts(self,dna_reads) :
        '''
        Input - list of dna reads
        Output - dictionary with key as read and value as the no. of times it occurs
        '''
        #start code here
        Read_counts = {}
        for i in range(len(dna_reads)):
            if dna_reads[i] not in Read_counts.keys():
                Read_counts[dna_reads[i]] = 1
            elif dna_reads[i] in Read_counts.keys():
                Read_counts[dna_reads[i]] = Read_counts[dna_reads[i]] + 1
        return Read_counts
        #end code here

    def parse_reads_pac(self,reads_pac):
        '''
        Input - pac bio reads file as a string
        Output - list of dna reads
        '''
        #start code here
        index = -1
        DNA_reads = []
        Splitlist_reads_pac = reads_pac.split('\n')
        for i in range(len(Splitlist_reads_pac)-1):
            if Splitlist_reads_pac[i][0] == ">":
                index = index + 1
                #print("i is",i)
                #print("index is", index)
                '''if index >= 1:
                    print("DNA_reads[index] is",DNA_reads[index-1])'''
                DNA_reads.append("")
            elif Splitlist_reads_pac[i][0] != ">":
                DNA_reads[index] = DNA_reads[index] + Splitlist_reads_pac[i]                
        return DNA_reads
        #end code here