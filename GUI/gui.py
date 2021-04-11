import collections
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import pylab
import math
import streamlit as st
import time
import os
from scipy.stats import spearmanr,kendalltau,pearsonr
from timer import Timer
import numpy as np
from matplotlib.backends.backend_agg import RendererAgg


st.set_option('deprecation.showPyplotGlobalUse', False)


class CGR():
    K = 0
    c = None
    h = ""
    Data = ""
    N=2
    def param(self,k,n):
        self.K = k
        self.N = n

    def read_fasta(self,loc):
        f = open(loc)
        s1 = f.read()
        data = "".join(s1.split("\n")[1:])
        head = "".join(s1.split("\n")[0:1])
        return data,head
    
    def count_kmers(self,sequence, k):
        d = collections.defaultdict(int)
        for i in range(len(self.Data)-(k-1)):
            d[sequence[i:i+k]] +=1
        d.pop("N",None)
        return d

    def probabilities(self,kmer_count, k):
        probabilities = collections.defaultdict(float)
        N = len(self.Data)
        for key, value in kmer_count.items():
            probabilities[key] = float(value) / (N - k + 1)
        return probabilities

    def chaos_game_representation(self,probabilities, k):
        array_size = int(math.sqrt(4**k))
        chaos = []
        for i in range(array_size):
            chaos.append([0]*array_size)
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
        for key, value in probabilities.items():
            for char in key:
                if char == "T":
                    posx +=  maxx/2
                elif char == "C":
                    posy += maxy/2
                elif char == "G":
                    posx += maxx/2
                    posy += maxy/2
                maxx /= 2
                maxy /= 2

            chaos[int(posy-1)][int(posx-1)] = value
            maxx = array_size
            maxy = array_size
            posx = 1
            posy = 1
        return chaos

    def load_fasta(self,loc):
        data,head = self.read_fasta(loc)
        self.Data = data
        t = Timer()
        t.start()
        f4 = self.count_kmers(data, self.K)
        f4_prob = self.probabilities(f4, self.K)
        chaos_k4 = self.chaos_game_representation(f4_prob, self.K)
        s = t.stop()
        self.c = chaos_k4
        self.h = head
        return chaos_k4,s
    
    def show(self):
        pylab.figure(figsize=(12,12))
        pylab.title('CGR of '+str(self.K)+'-mers for '+self.h[2:])
        pylab.imshow(self.c, cmap=cm.gray_r)
        ax = pylab.gca()
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        pylab.show()


def file_selector(folder_path='.'):
    filenames = os.listdir(folder_path)
    selected_filename = st.selectbox('Select a file', filenames)
    return os.path.join(folder_path, selected_filename)

def file_selector_1(folder_path='.'):
    filenames = os.listdir(folder_path)
    selected_filename = st.selectbox('Select comparision file', filenames)
    return os.path.join(folder_path, selected_filename)

if __name__ == '__main__':
        _lock = RendererAgg.lock
        st.title("Chaos Game Representation of Sequences")
        st.write(" ________ ")
        K = st.sidebar.slider(
        'Select K Value',
        2, 10, 6
        )
        A = CGR()
        B = CGR()

        A.param(K,2)
        B.param(K,2)
        #loc = st.file_uploader("Upload Fasta File",type=["fasta"])
        filename = file_selector()
        filename_1 = file_selector_1()

        cg,t1 = A.load_fasta(filename)
        cg_1,t2 = B.load_fasta(filename_1)
        st.write("_________")
        st.write(A.h)
        row3_space1, row3_1, row3_space2, row3_2, row3_space3 = st.beta_columns(
    (.1, 1, .1, 1, .1))
        with row3_1, _lock:
            A.show()
            st.pyplot()
        with row3_2, _lock:
            t = np.array(cg)
            plt.matshow(t)
            ax = plt.gca()
            ax.axes.xaxis.set_visible(False)
            ax.axes.yaxis.set_visible(False)
            plt.show()
            st.pyplot()
        st.write("_________")
        st.write(B.h)
        row4_space1, row4_1, row4_space2, row4_2, row4_space3 = st.beta_columns(
    (.1, 1, .1, 1, .1))
        
        with row4_1, _lock:
            B.show()
            st.pyplot()
        with row4_2, _lock:
            t_1 = np.array(cg)
            plt.matshow(t_1)
            ax = plt.gca()
            ax.axes.xaxis.set_visible(False)
            ax.axes.yaxis.set_visible(False)
            plt.show()
            st.pyplot()
        st.write("_________")

        
        a = max(cg[0])
        b = max(cg[1])
        CG = list(map(lambda x,y:((x*(1/a))+(y*(1/b))*10)**0.5, cg[0],cg[1]))

        a = max(cg_1[0])
        b = max(cg_1[1])
        CG_1 = list(map(lambda x,y:((x*(1/a))+(y*(1/b))*10)**0.5, cg_1[0],cg_1[1]))

        st.sidebar.text("Accuracy:")
        add_selectbox = st.sidebar.selectbox(
        'Similarity Method',
        ('Spearmans correlation', 'Kendalls tau','Pearsons correlation')
        )
        if(add_selectbox=='Spearmans correlation'): 
            corr, _ = spearmanr(CG,CG_1)
            st.sidebar.text(" ")
            st.sidebar.text(f'Spearmans correlation: %.3f' % corr)
        elif(add_selectbox=='Kendalls tau'):
            corr, _ = kendalltau(CG,CG_1)
            st.sidebar.text(" ")
            st.sidebar.text(f'Kendalls tau: %.3f' % corr)
        else:
            corr, _ = pearsonr(CG,CG_1)
            st.sidebar.text('Pearsons correlation: %.3f' % corr)  
        
        f = open(filename)
        s1 = f.read()
        data_1 = "".join(s1.split("\n")[1:])
        f = open(filename_1)
        s1 = f.read()
        data_2 = "".join(s1.split("\n")[1:])
        sum=0
        t=Timer()
        t.start()
        for i in range(len(data_2)):
            if data_1[i]==data_2[i]:
                sum+=1
        sum/=len(data_1)

        st.sidebar.text("for loop "+str(sum)+"\n")

        st.sidebar.text(".")
        st.sidebar.text("Run-Time:")
        st.sidebar.text("For loop-"+t.stop())
        st.sidebar.text("CGR - "+t1)
        st.sidebar.text(".")
        st.sidebar.text("Credits:")

        st.sidebar.text(""" 
            Team - O
            .    Aaditya Jain 
            .    Anirudh Bhaskar
            .    Ankam Srikanth
            .    Rohith Ramakrishnan
        """)
