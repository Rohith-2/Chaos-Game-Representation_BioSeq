# Chaos-Game-Representation_BioSeq
## AIE-19BIO211
### Intelligence In Biological Systems - 4
Representation of Bio Sequences via Chaos Game and using the same to find similarities

### Data
Each excel file is a collection of all sequences in the respective family. 

### GUI
https://share.streamlit.io/rohith-2/chaos-game-representation_bioseq/stream.py  

<img width="1327" alt="Screenshot 2021-05-01 at 11 13 17 AM" src="https://user-images.githubusercontent.com/55501708/116772668-483cca00-aa6e-11eb-8727-06e945e1ea9a.png">


### Gene-Similarity
CGR Matrix is a 2D matrix => (x,y) which consists of normalised value ranging from 0 to 1, which depicts the intensity of a color at any given (x,y)  
The first two rows are considered for similarity measurement:  
```Python
a = max(first row of cgr matrix of SEQ_1)
b = max(second row of cgr matrix of SEQ_1)
CG = list(map(lambda x,y:((x*(1/a))+(y*(1/b))*10)**0.5, cg[0],cg[1]))
a = max(first row of cgr matrix of SEQ_2)
b = max(second row of cgr matrix of SEQ_2)
CG_1 = list(map(lambda x,y:((x*(1/a))+(y*(1/b))*10)**0.5, cg_1[0],cg_1[1]))
```

CG and CG_1 will be vectors which can be utilised for measuring similarity via Spearmans correlation. 

### Running the GUI
The GUI is universally accesible via the above mentioned link to run in locally :  
```
git clone https://github.com/Rohith-2/Chaos-Game-Representation_BioSeq.git
pip install -r requirements.txt  
cd Chaos-Game-Representation_BioSeq/GUI/
streamlit run gui_v2.py 
```
