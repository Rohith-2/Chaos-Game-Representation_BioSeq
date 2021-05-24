# Chaos-Game-Representation_BioSeq  

### Representation of Bio Sequences via Chaos Game and using the same to find similarities.   
__Medium : https://rrohith2001.medium.com/chaos-game-representation-of-genetic-sequences-e0e6bdcfaf6c__  
<hr style=\"border:0.5px solid gray\"> </hr>  

### Data    
Each excel file is a collection of all sequences in the respective family. 

<hr style=\"border:0.5px solid gray\"> </hr>  

### GUI  
https://share.streamlit.io/rohith-2/chaos-game-representation_bioseq/stream.py  

![Screenshot 2021-05-18 at 7 41 30 PM](https://user-images.githubusercontent.com/55501708/118666828-19b24380-b811-11eb-8ad6-ac894089cf01.png)

<hr style=\"border:0.5px solid gray\"> </hr>  
  
### Gene-Similarity  
CGR Matrix is a 2D matrix => (x,y) which consists of normalised value ranging from 0 to 1, which depicts the intensity of a color at any given (x,y)  
The first two rows are considered for similarity measurement:  
```Python
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
```

cgr_vec and cgr_vec_2 will be vectors which can be utilised for measuring similarity via Spearmans correlation. 
<hr style=\"border:0.5px solid gray\"> </hr>  

### Running the GUI  
The GUI is universally accesible via the above mentioned link to run in locally :  
```
git clone https://github.com/Rohith-2/Chaos-Game-Representation_BioSeq.git
pip install -r requirements.txt  
cd Chaos-Game-Representation_BioSeq/GUI/
streamlit run gui_v2.py 
```
