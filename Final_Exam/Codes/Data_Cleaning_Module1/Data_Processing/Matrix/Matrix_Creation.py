#Matrix creation
import numpy as np
import json
import itertools
from guess_language import guess_language
import re
import scipy.io


def is_ascii(s):
    return all(ord(c) < 128 for c in s)

filename = "JD.txt"
target = open(filename, 'r')

file_contents = target.read()

f = file_contents.split()

index = []
for i in range(len(file_contents)):
     index.append(i)

dictionary = dict(zip(f,index))

with open('Java Developer.json') as json_data:
    d = json.load(json_data)
    
    no_users = len(d)
    M = np.zeros([len(dictionary),no_users])
    
    for i in range(len(d)):
    		D = d[i]['Skills']
    	
		G = D.split(",")
		H = "php"
		H1= 'lua'
		H2= 'jsf'
		H3= 'css'
		H4= 'jquery'
		H5 = 'gwt'
		H6 = 'musql'
		H7 = 'mongodb'
		H8 = 'jira'
		H9 = 'angularjs' 
		H10 = 'junit'
		H11 = 'seleniun'
		H12 = 'akka-http'
		H13 = 'soa'
		H14 = 'java/j2ee'
		H15 = 'Servlets'
		H16 = 'Spring MVC'
		
		for j in range(len(G)):
			if(is_ascii(G[j])):
				K = re.sub('\(.*?\)','', str(G[j]))
			
				if(''.join(K.split()).lower()=='java'):
						M[0][i]=1
				elif (''.join(K.split()).lower()==H):
						M[dictionary[H]][i]=1
			   	elif (''.join(K.split()).lower()==H1):
						M[dictionary[H1]][i]=1 	
				elif (''.join(K.split()).lower()==H2):
						M[dictionary[H2]][i]=1
				elif (''.join(K.split()).lower()==H3):
						M[dictionary[H3]][i]=1
				elif (''.join(K.split()).lower()==H4):
						M[dictionary[H4]][i]=1
				elif (''.join(K.split()).lower()==H5):
						M[dictionary[H5]][i]=1
				elif (''.join(K.split()).lower()==H6):
						M[dictionary[H6]][i]=1
				elif (''.join(K.split()).lower()==H7):
						M[dictionary[H7]][i]=1
				elif (''.join(K.split()).lower()==H8):
						M[dictionary[H8]][i]=1
				elif (''.join(K.split()).lower()==H9):
						M[dictionary[H9]][i]=1 	
				elif (''.join(K.split()).lower()==H10):
						M[dictionary[H10]][i]=1
				elif (''.join(K.split()).lower()==H11):
						M[dictionary[H11]][i]=1		 	
				elif (guess_language(''.join(K.split()).upper())=='en'):
						M[dictionary[''.join(K.split()).lower()]][i]=1				
				elif (guess_language(''.join(K.split()).lower())=='en'):
						M[dictionary[''.join(K.split()).lower()]][i]=1
							
				elif (guess_language(K)=='en'):
						M[dictionary[''.join(K.split()).lower()]][i]=1

scipy.io.savemat('output.mat', mdict={'M': M})


