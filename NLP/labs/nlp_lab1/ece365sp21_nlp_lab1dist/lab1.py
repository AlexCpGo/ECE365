from nltk.stem import WordNetLemmatizer
from nltk.stem import PorterStemmer
import matplotlib.pyplot as plt
from nltk.corpus import udhr
import re

def get_freqs(corpus, puncts):
    freqs = {}
    ### BEGIN SOLUTION
    for i in puncts:
        if i in corpus:
            corpus=corpus.replace(i,' ')
    corpus=re.sub('[0-9]', ' ', corpus)
    corpus=corpus.lower()
    data=corpus.split( )
    ### END SOLUTION
    for i in data:
        if i not in freqs:
            freqs[i]=1
        else:
            freqs[i]+=1
    return freqs

def get_top_10(freqs):
    top_10 = []
    ### BEGIN SOLUTION
    SortedFreqs=sorted(freqs.items(), key=lambda item: item[1], reverse=True)
    top_10_freqs = SortedFreqs[:10]
    top_10=[i[0] for i in top_10_freqs]
    ### END SOLUTION
    return top_10

def get_bottom_10(freqs):
    bottom_10 = []
    ### BEGIN SOLUTION   
    SortedFreqs=sorted(freqs.items(), key=lambda item: item[1])
    bottom_10_freqs=SortedFreqs[:10]
    bottom_10=[i[0] for i in bottom_10_freqs]
    return bottom_10
    ### END SOLUTION
    return bottom_10

def get_percentage_singletons(freqs):
    ### BEGIN SOLUTION
    ### END SOLUTION
    pass

def get_freqs_stemming(corpus, puncts):
    ### BEGIN SOLUTION
    ### END SOLUTION
    pass

def get_freqs_lemmatized(corpus, puncts):
    ### BEGIN SOLUTION
    ### END SOLUTION
    pass

def size_of_raw_corpus(freqs):
    ### BEGIN SOLUTION
    ### END SOLUTION
    pass

def size_of_stemmed_raw_corpus(freqs_stemming):
    ### BEGIN SOLUTION
    ### END SOLUTION
    pass

def size_of_lemmatized_raw_corpus(freqs_lemmatized):
    ### BEGIN SOLUTION
    ### END SOLUTION
    pass

def percentage_of_unseen_vocab(a, b, length_i):
    ### BEGIN SOLUTION
    ### END SOLUTION
    pass

def frac_80_perc(freqs):
    ### BEGIN SOLUTION
    ### END SOLUTION
    pass

def plot_zipf(freqs):
    ### BEGIN SOLUTION
    ### END SOLUTION
    plt.show()  # put this line at the end to display the figure.

def get_TTRs(languages):
    TTRs = {}
    for lang in languages:
        words = udhr.words(lang)
        ### BEGIN SOLUTION
        ### END SOLUTION
    return TTRs

def plot_TTRs(TTRs):
    ### BEGIN SOLUTION
    ### END SOLUTION
    plt.show()  # put this line at the end to display the figure.
