import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
class Lab4(object):
    
    def expectation_maximization(self,read_mapping,tr_lengths,n_iterations) :
        #start code here
        H=np.zeros((len(tr_lengths),1+n_iterations))
        z=np.zeros((len(read_mapping),len(tr_lengths)))
        for i in range(len(tr_lengths)):
            H[i][0]=1/len(tr_lengths)
        for j in range(n_iterations):
            for i in range(len(read_mapping)):
                for k in range(len(tr_lengths)):
                    if k in read_mapping[i]:
                        s=0
                        for e in read_mapping[i]:
                            s+=H[e][j]
                        z[i][k]=H[k][j]/s
                    else:
                        z[i][k]=0
            the=np.zeros(len(tr_lengths))
            for k in range(len(tr_lengths)):
                the[k]=z.sum(axis=0)[k]/len(read_mapping)
            for k in range(len(tr_lengths)):
                H[k][j+1]=the[k]/tr_lengths[k]/(the/tr_lengths).sum()
        return H
                
        #end code here

    def prepare_data(self,lines_genes) :
        '''
        Input - list of strings where each string corresponds to expression levels of a gene across 3005 cells
        Output - gene expression dataframe
        '''
        #start code here
        GeneDict={}
        gene_idx=0
        for line in lines_genes:
            gene_k_array=np.round(np.log(np.array(line.split(' ')).astype("float")+1),5)
            GeneDict["Gene_"+str(gene_idx)]=gene_k_array
            gene_idx=gene_idx+1
        df=pd.DataFrame(GeneDict)
        return df
        
        #end code here
    
    def identify_less_expressive_genes(self,df) :
        '''
        Input - gene expression dataframe
        Output - list of column names which are expressed in less than 25 cells
        '''
        #start code here
        GeneLis=[]
        col=df.columns
        mat=df.values>0
        s=mat.sum(axis=0)
        for i in range(len(s)):
            if s[i]<25:
                GeneLis.append(col[i])
        return GeneLis        
                
        #end code here
    
    
    def perform_pca(self,df) :
        '''
        Input - df_new
        Output - numpy array containing the top 50 principal components of the data.
        '''
        #start code here
        pca=PCA(n_components=50,random_state=365)
        result=np.round(pca.fit_transform(df.values),5)
        return result
        #end code here
    
    def perform_tsne(self,pca_data) :
        '''
        Input - pca_data
        Output - numpy array containing the top 2 tsne components of the data.
        '''
        #start code here
        result=TSNE(n_components=2,perplexity=50, random_state=1000).fit_transform(pca_data)
        return result
        #end code here