import collections
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import pylab
import math
import pandas as pd
import streamlit as st
from scipy.stats import spearmanr,kendalltau,pearsonr
import time
import numpy as np
import os
from matplotlib.backends.backend_agg import RendererAgg
import cv2
from math import log10, sqrt
from SSIM_PIL import compare_ssim
from PIL import Image
from skimage.metrics import structural_similarity as ssim
import matplotlib 
matplotlib.use("Agg") 
st.set_option('deprecation.showPyplotGlobalUse', False)



l = os.getcwd()
if(l!='/app/chaos-game-representation_bioseq/data/'):
     os.chdir('/app/chaos-game-representation_bioseq/data/')

class TimerError(Exception):
     """A custom exception used to report errors in use of Timer class"""
 
class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None
        return (f"Elapsed time: {elapsed_time:0.4f} seconds")


def PSNR(original, compressed):
    mse = np.mean((original - compressed) ** 2)
    if(mse == 0):  # MSE is zero means no noise is present in the signal .
                  # Therefore PSNR have no importance.
        return 100
    max_pixel = 255.0
    psnr = 20 * log10(max_pixel / sqrt(mse))
    return psnr

class CGR():
    K = 0
    c = None
    h = ""
    Data = ""
    N=2
    F=""
    i = 0
    def __init__(self,a):
        self.i=a

    def param(self,k,n):
        self.K = k
        self.N = n

    def get_family(self):
        return self.F

    def read_fasta(self,family,locz):
        x = pd.read_excel(family,engine = 'openpyxl',header=None)
        x.set_index(0, inplace = True)
        data = x.loc[locz].tolist()
        self.F = family
        return data,locz
    
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

    def load_fasta(self,data,head):
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
        pylab.imshow(self.c, cmap=cm.gray_r)#,filterrad=10,vmin=0,vmax=2)
        ax = pylab.gca()
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        pylab.axis("off")
        pylab.savefig(str(self.i)+".PNG",bbox_inches='tight', pad_inches=0)
        pylab.show()


def file_selector(F):
    x = pd.read_excel(F ,engine = 'openpyxl',header=None)
    acc = x[0].values.tolist()
    def format_func(option):
        return option.split(',')[0]
    selected_filename = st.selectbox('Select a Sequence', acc,format_func=format_func)
    print(selected_filename)
    x.set_index(0, inplace = True)
    a = x.loc[selected_filename].tolist()
    return a[0],selected_filename

def file_selector_1(F):
    x = pd.read_excel(F ,engine = 'openpyxl',header=None)
    acc = x[0].values.tolist()
    def format_func(option):
        return option.split(',')[0]
    selected_filenam = st.selectbox('Select comparision Sequence', acc,format_func=format_func)
    x.set_index(0, inplace = True)
    a = x.loc[selected_filenam].tolist()
    return a[0],selected_filenam

if __name__ == '__main__':
        _lock = RendererAgg.lock
        st.title("Chaos Game Representation of Sequences")
        st.write(" ________ ")
        K = st.sidebar.slider(
        'Select K Value',
        4, 12, 6
        )
        A = CGR("F")
        B = CGR("S")

        A.param(K,2)
        B.param(K,2)
        def format_func(option):
            return option[:-5]

        filenames = os.listdir()
        family = st.selectbox('Select Family', filenames,format_func=format_func)
        filename,h1 = file_selector(family)

        family1 = st.selectbox('Select Comparision Family', filenames,format_func=format_func)
        filename_1,h2 = file_selector_1(family1)
        

        cg,t1 = A.load_fasta(filename,h1)
        cg_1,t2 = B.load_fasta(filename_1,h2)
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

        t = Timer()
        t.start()
        imageA = cv2.imread("F.PNG")
        imageB = cv2.imread("S.PNG")
        grayA = cv2.cvtColor(imageA, cv2.COLOR_BGR2GRAY)
        grayB = cv2.cvtColor(imageB, cv2.COLOR_BGR2GRAY)
        (score,_) = ssim(grayA, grayB, full=True)
        ssim = t.stop()

        value = PSNR(grayA, grayB)
        st.write("About PSNR:")
        st.write("PSNR is most commonly used to estimate the efficiency of compressors, filters, etc. The larger the value of PSNR, the more efficient is a corresponding compression or filter method.")
        st.write("_________")
        st.write("About SSIM:")
        st.write("The Structural Similarity Index (SSIM) is a perceptual metric that quantifies image quality degradation* caused by processing such as data compression or by losses in data transmission.")
        st.write("_________")
        st.write("Packages used for SSIM:")
        st.write("SSIM-PIL : https://pypi.org/project/SSIM-PIL/ & https://scikit-image.org/docs/dev/auto_examples/transform/plot_ssim.html")
        st.write("_________")
        st.write("Algorithm for correlation:")
        st.markdown('''
            ```python
          cgr_vec = Empty Vector()
          for i <- cgr matrix of SEQ_1 # i iterates row wise
            a = max(i)
            new_row = i/a # Element-Wise Division
            cgr_vector = cgr_vector + new_row

          cgr_vec_2 = Empty Vector()
          for i <- cgr matrix of SEQ_2 # i iterates row wise
            a = max(i)
            new_row = i/a # Element-Wise Division
            cgr_vector_2 = cgr_vector_2 + new_row

          Correlation(cgr_vec,cgr_vec_2)
            ''')

        t = Timer()
        t.start()
        image1 = Image.open('F.PNG')
        image2 = Image.open('S.PNG')
        value1 = compare_ssim(image1, image2)
        pil = t.stop()

        t = Timer()
        t.start()
        cg = np.array(cg)
        val_1=np.zeros(cg.shape[1])
        for i in cg:
            if max(i)!= 0: 
                a = max(i) 
            else:
                a = 0.0001
            x = i*(1/a)
            val_1 += x
        CG = val_1**0.5

        cg_1 = np.array(cg_1)
        val_2=np.zeros(cg_1.shape[1])
        for i in cg_1:
            if max(i)!= 0: 
                a = max(i) 
            else:
                a = 0.0001
            x = i*(1/a)
            val_2 += x
        CG_1 = val_2**0.5
        com = t.stop()

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
    
        st.sidebar.text("SSIM Score: "+str(score))
        st.sidebar.text("SSIM-PIL value: "+str(value1))
        st.sidebar.text("PSNR value: "+str(value)+"dB")
        st.sidebar.text(".")
        st.sidebar.text("Run-Time:")
        st.sidebar.text("CGR - "+str(t1))
        st.sidebar.text("Correlation - "+str(com))
        st.sidebar.text("SSIM - "+str(ssim))
        st.sidebar.text("SSIM-PIL - "+str(pil))
        st.sidebar.text(".")
        st.sidebar.text("Credits:")

        st.sidebar.text(""" 
            Team - O
            .    Aaditya Jain 
            .    Anirudh Bhaskar
            .    Ankam Srikanth
            .    Rohith Ramakrishnan
        """)
        st.write("_________")
        st.text("Meduim Article:")
        st.write("https://rrohith2001.medium.com/chaos-game-representation-of-genetic-sequences-e0e6bdcfaf6c")
