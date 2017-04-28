#Data Cleaning 

#Here issues happens when skills like lua is not recognize in english dictionary 
#Skills with more than 1 year experince will be considerd 
import json
from guess_language import guess_language
import re
from sys import argv


def is_ascii(s):
    return all(ord(c) < 128 for c in s)


filename = "Software Developer - Backend Companies.txt"
target = open(filename, 'w')

with open('Software Developer - Backend.json') as json_data:
    d = json.load(json_data)
    
    for i in range(len(d)):
    	D = d[i]['Work-Experience']['Company']
	if(is_ascii(D)):
		target.write(D)
		
	target.write("\n")

	
