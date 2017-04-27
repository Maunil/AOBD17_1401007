#Skills Ordering 
import numpy as np
import pickle
import itertools
import csv
import operator 

def unique_file(input_filename, output_filename):
    input_file = open(input_filename, 'r')
    file_contents = input_file.read()
    input_file.close()
    duplicates = []
    word_list = file_contents.split()
    file = open(output_filename, 'w')
    
    for word in word_list:
        if word not in duplicates:
             duplicates.append(word)    
             #file.write(word)
             #file.write("\n")

    index = []
    for i in range(len(duplicates)):
        index.append(1)
    #index = np.zeros([len(duplicates)])
    dictionary = dict(itertools.izip(duplicates,index))  
    
    for word in word_list:
        if word  in duplicates:
             #dictionary.get(word)                
             dictionary[word] = dictionary[word]+1   
     
    sorted_x = sorted(dictionary.items(), key=operator.itemgetter(1))
    
    for i in range(len(dictionary)):
           file.write(sorted_x[len(dictionary)-1-i][0])
           file.write("\n")
    
    file.close()

unique_file('Java Developer.txt','JD.txt')

