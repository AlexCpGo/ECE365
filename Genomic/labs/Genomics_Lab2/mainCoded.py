import numpy as np
from collections import OrderedDict

class Lab2(object):
    
    def smith_waterman_alignment(self,s1,s2,penalties) :
        '''
        Input - two sequences and a dictionary with penalities for match, mismatch and gap
        Output - an integer value which is the maximum smith waterman alignment score
        '''
        #start code here
        # start_time1 = time.time()
        # H store the score table of the smith waterman algorithm with the size of L(s1) * L(s2) 
        H = np.zeros((len(s1)+1,len(s2)+1))
        # Fill the first column with zeros
        # H.append(np.zeros(len(s2)+1))
        match = penalties['match']
        mismatch = penalties['mismatch']
        gap = penalties['gap']
        MaxScore = 0
         
        # s1->row; s2->column
        for i in range(1,len(s1)+1): 
            H[i][0] = 0
            for j in range(1,len(s2)+1):
                if s1[i-1] == s2[j-1]:
                    H_change = H[i-1][j-1]+match
                elif s1[i-1] != s2[j-1]:
                    H_change = H[i-1][j-1]+mismatch

                H[i][j] = max(H_change, H[i-1][j]+gap, H[i][j-1]+gap, 0)
                if H[i][j] > MaxScore:
                    MaxScore = H[i][j]
        #print("--- %s seconds ---" % (time.time() - start_time1))
        return MaxScore
        #end code here

    def print_smith_waterman_alignment(self,s1,s2,penalties) :
        '''
        Input - two sequences and a dictionary with penalities for match, mismatch and gap
        Output - a tuple with two strings showing the two sequences with '-' representing the gaps
        '''
        #start code here
        #start_time2 = time.time()
        match = penalties['match']
        mismatch = penalties['mismatch']
        gap = penalties['gap']

        len_x=len(s1)
        len_y=len(s2)      
        H = np.zeros((len_x+1,len_y+1))        
        MaxScore=0

        for i in range(len_x+1):
            H[i][0]=0
        for j in range(len_y+1):
            H[0][j]=0
        for i in range(len_x):
            for j in range(len_y):
                if s1[i]==s2[j]:
                    H_change = H[i][j] + match
                else:
                    H_change = H[i][j] + mismatch
                Score = max(H[i][j+1] + gap, H[i+1][j] + gap, H_change,0)
                
                if Score > MaxScore:
                    MaxScore = Score
                    Max_i = i+1
                    Max_j = j+1
                H[i+1][j+1] = Score

        alignment = ['','']     
        i = Max_i
        j = Max_j
        while H[i][j]:
            if s1[i-1] == s2[j-1]:
                alignment[0] = s1[i-1] + alignment[0]
                alignment[1] = s2[j-1] + alignment[1]
                i = i - 1
                j = j - 1
            elif H[i][j-1] + gap == H[i][j]:
                alignment[0] = "-" + alignment[0]
                alignment[1] = s2[j-1] + alignment[1]
                j = j - 1
            elif H[i-1][j] + gap == H[i][j]:
                alignment[0] = s1[i-1] + alignment[0]
                alignment[1] = "-" + alignment[1]
                i = i - 1
            else:
                alignment[0] = s1[i-1] + alignment[0]
                alignment[1] = s2[j-1] + alignment[1]
                i = i - 1
                j = j - 1
        #print("--- %s seconds ---" % (time.time() - start_time2))      
        return tuple(alignment)
        #end code here

    def find_exact_matches(self,list_of_reads,genome):
       
        '''
        Input - list of reads of the same length and a genome fasta file (converted into a single string)
        Output - a list with the same length as list_of_reads, where the ith element is a list of all locations (starting positions) in the genome where the ith read appears. The starting positions should be specified using the "chr2:120000" format
        '''
        #start code here
        #start_time3 = time.time()
        NewGenome= genome.split('>')[1:]
        GenomeData=[]
        GenomeDict={}
        N = []
        for i in range(0,10):
            N.append(str(i))

        # Pre_processing with the original read data
        for rd in NewGenome:
            rd=rd.replace('\n','')
            for l in N:
                rd=rd.replace(l,'')
            GenomeData.append(rd.replace('chr',''))
        
        # Build the dictionary of genome
        length=len(list_of_reads[0])
        for i in range(len(GenomeData)):
            if len(GenomeData[i])>=length:
                for j in range(len(GenomeData[i])-length+1):
                    if GenomeData[i][j:j+length] not in GenomeDict:
                        label='chr'+str(i+1)+':'+str(j+1)
                        GenomeDict[GenomeData[i][j:j+length]]=[label]
                    else:
                        label='chr'+str(i+1)+':'+str(j+1)
                        GenomeDict[GenomeData[i][j:j+length]].append(label)

        # Match the list_of_read with Genome Dictionary
        MatchList=[]
        for m in list_of_reads:
            if m in GenomeDict:
                MatchList.append(GenomeDict[m])
            else:
                MatchList.append([])
        #print("--- %s seconds ---" % (time.time() - start_time3))       
        return MatchList
        
        #end code here
       
    
    def find_approximate_matches(self,list_of_reads,genome):
        '''
        Input - list of reads of the same length and a genome fasta file (converted into a single string)
        Output -  a list with the same length as list_of_reads, where the ith element is a list of all locations (starting positions) in the genome which have the highest smith waterman alignment score with ith read in list_of_reads
        '''
        #start code 
        #start_time4 = time.time()
        N = []
        for i in range(0,10):
            N.append(str(i))
        New_genome= genome.split('>')[1:]
        GenomeData=[]
        GenomeDict={}
        Matchlist=[]
        MaxScore=0
        max_list=[]

        # Pre_processing on the read genome data
        penalties={'match':1,'mismatch':-1,'gap':-1}
        for rd in New_genome:
            rd=rd.replace('\n','')
            for e in N:
                rd=rd.replace(e,'')
            GenomeData.append(rd.replace('chr',''))
        
        length=int(len(list_of_reads[0])/4)
 
        # Build up the Genome Dictionary
        for i in range(len(GenomeData)):
            if len(GenomeData[i])>=length:
                for j in range(len(GenomeData[i])-length+1):
                    if GenomeData[i][j:j+length] not in GenomeDict:
                        label='chr'+str(i+1)+':'+str(j+1)
                        GenomeDict[GenomeData[i][j:j+length]]=[label]
                    else:
                        label='chr'+str(i+1)+':'+str(j+1)
                        GenomeDict[GenomeData[i][j:j+length]].append(label)

        map=[[0]*(len(list_of_reads[0])-length+1) for _ in range(len(list_of_reads))]
        for i in range(len(list_of_reads)):
            for j in range(len(list_of_reads[0])-length+1):
                segment=list_of_reads[i][j:int(j+len(list_of_reads[0])/4)]
                if segment in GenomeDict:
                    map[i][j]=GenomeDict[segment]

        # Use the read data to match with the dictionary
        for i in range(len(list_of_reads)):
            max_list=[]
            for j in range(len(map[i])):
                if map[i][j]!=0:
                    for k in range(len(map[i][j])):
                        msg=map[i][j][k]
                        channel=int(msg.split(':')[0].split('chr')[1])
                        idx=int(msg.split(':')[-1])
                        x=list_of_reads[i]
                        if idx-j-1+len(x)<=len(GenomeData[channel-1]):
                            y=GenomeData[channel-1][idx-1-j:idx-j-1+len(x)]

                            score=self.smith_waterman_alignment(x,y,penalties)
                        
                            if len(max_list)==0:
                                MaxScore=score
                                if ('chr'+str(channel)+':'+str(idx-j)) not in max_list:
                                    max_list.append('chr'+str(channel)+':'+str(idx-j))
                            else:
                                if len(max_list)>0 and score==MaxScore:
                                    if ('chr'+str(channel)+':'+str(idx-j)) not in max_list:
                                        max_list.append('chr'+str(channel)+':'+str(idx-j))
                                elif score>MaxScore:
                                    MaxScore=score
                                    max_list=[]
                                    if ('chr'+str(channel)+':'+str(idx-j)) not in max_list:
                                        max_list.append('chr'+str(channel)+':'+str(idx-j))

            Matchlist.append(max_list)
        #print("--- %s seconds ---" % (time.time() - start_time4))   
        return Matchlist
        #end code here