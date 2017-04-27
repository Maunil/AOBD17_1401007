#Data Cleaning 
 
import json
from guess_language import guess_language
import re
from sys import argv


def is_ascii(s):
    return all(ord(c) < 128 for c in s)


filename = "Junior Software Developer Companies.txt"
target = open(filename, 'w')

with open('Junior Software Engineer.json') as json_data:
    d = json.load(json_data)
    
    for i in range(len(d)):
    	D = d[i]['Work-Experience']['Company']
	if(is_ascii(D)):
		target.write(D)
		
	target.write("\n")

	
