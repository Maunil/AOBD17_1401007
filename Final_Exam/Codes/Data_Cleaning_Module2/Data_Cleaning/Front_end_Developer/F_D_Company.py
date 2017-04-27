#Data Cleaning 
#Company Data Extraction 

import json
from guess_language import guess_language
import re
from sys import argv


def is_ascii(s):
    return all(ord(c) < 128 for c in s)


filename = "Front End Developer Companies.txt"
target = open(filename, 'w')

with open('Front End Developer.json') as json_data:
    d = json.load(json_data)
    
    for i in range(len(d)):
    	D = d[i]['Work-Experience']['Company']
    	print D
	if(is_ascii(D)):
		target.write(''.join(D.split()).lower())
		target.write("\n")
		
	target.write("\n\n")

	
