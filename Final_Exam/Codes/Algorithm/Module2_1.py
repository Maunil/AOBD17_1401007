import re
import os
import json
import pickle

def write(data, outfile):
    with open(outfile, 'w') as fp:
        json.dump(data, fp, indent=4)

def read(filename):
    with open(filename, 'r') as fp:
        data = json.load(fp)
    return data

def year_edit(year_dict,year,skill,prof):

    if prof in year_dict.keys():
        year_dict[prof][skill] = year
    else:
        year_dict[prof] = {}
        year_dict[prof][skill] = year

    return year_dict


def company_processing(filepath,country,prof):
    path = ''.join([filepath,'/companies.json'])

    if os.path.exists(path):
        companies = read(path)
        if prof in list(companies.keys()):
            if country in list(companies[prof].keys()):
                return companies[prof][country]
            else:
                dict = {}
                dict[0] = "no jobs in "+country
                return dict
        else:
            dict = {}
            dict[0] = "no jobs in "+country
            return dict
    else:
        companies = {}
        for filename in os.listdir(filepath):
            c = 0
            string = ''
            with open(''.join([filepath,'/',filename])) as f:
                for line in f:
                    if c==0:
                        tmp = re.sub('\n','', line)
                        key = tmp
                    elif c>=2:
                        tmp = re.sub('\n','', line)
                        count = tmp.rsplit(':', 1)[0]
                        string = tmp.rsplit(':', 1)[1]
                        count = re.sub(' ','', count)

                        if key in companies.keys():
                            if count in companies[key].keys():
                                companies[key][count].append(string)
                            else:
                                companies[key][count] = [string]
                        else:
                            companies[key] = {}
                            companies[key][count] = [string]
                    c=c+1

        write(companies,path)
        if prof in list(companies.keys()):
            if country in list(companies[prof].keys()):
                return companies[prof][country]
            else:
                dict = {}
                dict[0] = "no jobs in "+country
                return dict
        else:
            dict = {}
            dict[0] = "no jobs in "+country
            return dict


folder = "My_Files_2"
json_folder = "JSON_Files_2"
new_dict = {} # Stores skill wise carriers
file_dict = {} # Stores carrier wise required skills
year_dict = {}
filepath_companies = "Companies"
additional_info_path = "Additional_Info"

new_dict_path = ''.join([json_folder,'/new_dict.json'])
file_dict_path = ''.join([json_folder,'/file_dict.json'])
year_dict_path = ''.join([json_folder,'/year_dict.json'])

if os.path.exists(new_dict_path):

    new_dict = read(new_dict_path)
    file_dict = read(file_dict_path)
    year_dict = read(year_dict_path)
else:
    for filename in os.listdir(folder):

        c = 0
        string = ''

        with open(''.join([folder,'/',filename])) as f:
            for line in f:
                if c==0:
                    string = re.sub('\n','', line)

                    if string not in file_dict.keys():
                        file_dict[string] = []
                elif c>=2:
                    tmp = re.sub('\n', '', line)
                    print tmp
                    vec = tmp.rsplit(':', 1)
                    year = vec[1]
                    tmp = vec[0]
                    year_dict = year_edit(year_dict,year,tmp,string)
                    if tmp in new_dict.keys():
                        new_dict[tmp].append(string)
                    else:
                        new_dict[tmp] = [string]

                    file_dict[string].append(tmp)

                c = c+1

    write(new_dict,new_dict_path)
    write(file_dict,file_dict_path)
    write(year_dict,year_dict_path)

# print file_dict

k = 0
skill_dict={}   # Stores Skills of user
with open('new_user_2.txt') as f:
    for line in f:
        if k==0:
            user_name = re.sub('\n','', line)
            if user_name not in skill_dict.keys():
                skill_dict[user_name] = []
        elif k==1:
            country = re.sub('\n','', line)
        elif k==2:
            carrier = re.sub('\n','', line)
        elif k>=4:
            tmp = re.sub('\n','', line)
            skill_dict[user_name].append(tmp)

        k = k+1

# print carrier_dict
# print skill_dict
username = user_name
os.system("clear")
num = []

li = []
cnt = 0
for item in file_dict[carrier]:
    if item not in skill_dict[username]:
        li.append(item)

print "Hello",username,":-)"
print "If you want to become", carrier, "you have to acquire following skills:\n"

for k in range(len(li)):
    print "    ->",li[k],", Min. Experience :",year_dict[carrier][li[k]][0]
print "\n****Who will hire you?****"
print "\nThe Companies that could hire you after becoming",carrier,"in",country,"are:\n"
dct = company_processing(filepath_companies,country,carrier)
for k in range(len(dct)):
    print "    ->",dct[k]
print "\n****Additional skills that may be helpful to you****"
print "->People working as",carrier,"have following additional skills:-\n"
tmp_path = ''.join([additional_info_path,'/',carrier,'.txt'])
if os.path.exists(tmp_path):
    with open(tmp_path) as f:
        for line in f:
            print line
