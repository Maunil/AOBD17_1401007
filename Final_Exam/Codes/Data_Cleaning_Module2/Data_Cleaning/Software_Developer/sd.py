#Data Cleaning 

#Here issues happens when skills like lua is not recognize in english dictionary 
#Skills with more than 1 year experince will be considerd 
import json
from guess_language import guess_language
import re
from sys import argv


def is_ascii(s):
    return all(ord(c) < 128 for c in s)


filename = "Software Developer.txt"
target = open(filename, 'w')

with open('Software Developer.json') as json_data:
    d = json.load(json_data)
    
    for i in range(len(d)):
    	D = d[i]['Skills']
    	print D

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
		#print G[j]
		if(is_ascii(G[j])):
			K = re.sub('\(.*?\)','', str(G[j]))
			
			#K = ''.join(K.split()).lower()
			if (''.join(K.split()).lower()==H):
					target.write(H)
					target.write("\n")
		   	elif (''.join(K.split()).lower()==H1):
					target.write(H1)
					target.write("\n") 	
			elif (''.join(K.split()).lower()==H2):
					target.write(H2)
					target.write("\n")
			elif (''.join(K.split()).lower()==H3):
					target.write(H3)
					target.write("\n")
			elif (''.join(K.split()).lower()==H4):
					target.write(H4)
					target.write("\n")
			elif (''.join(K.split()).lower()==H5):
					target.write(H5)
					target.write("\n")
			elif (''.join(K.split()).lower()==H6):
					target.write(H6)
					target.write("\n")
			elif (''.join(K.split()).lower()==H7):
					target.write(H7)
					target.write("\n")
			elif (''.join(K.split()).lower()==H8):
					target.write(H8)
					target.write("\n")
			elif (''.join(K.split()).lower()==H9):
					target.write(H9)
					target.write("\n") 	
			elif (''.join(K.split()).lower()==H10):
					target.write(H10)
					target.write("\n")
			elif (''.join(K.split()).lower()==H11):
					target.write(H11)
					target.write("\n")		 	
			elif (''.join(K.split()).lower()==H12):
					target.write(H12)
					target.write("\n")		 	
			
			elif (''.join(K.split()).lower()==H13):
					target.write(H13)
					target.write("\n")		 	
			elif (''.join(K.split()).lower()==H14):
					target.write(H14)
					target.write("\n")	
			elif (''.join(K.split()).lower()==H15):
					target.write(H15)
					target.write("\n")	
			elif (''.join(K.split()).lower()==H16):
					target.write(H16)
					target.write("\n")	
			
						 	
			elif (guess_language(''.join(K.split()).upper())=='en'):
				target.write(''.join(K.split()).lower())
				target.write("\n")
			elif (guess_language(''.join(K.split()).lower())=='en'):
				target.write(''.join(K.split()).lower())
				target.write("\n")	
			elif (guess_language(K)=='en'):
				target.write(''.join(K.split()).lower())
				target.write("\n")	

	target.write("\n\n")

