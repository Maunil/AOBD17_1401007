import re
import os
import json

def write(data, outfile):
    with open(outfile, 'w') as fp:
        json.dump(data, fp, indent=4)

def read(filename):
    with open(filename, 'r') as fp:
        data = json.load(fp)
    return data

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
            

folder = "My_Files"
json_folder = "JSON_Files"
new_dict = {} # Stores skill wise carriers
file_dict = {} # Stores carrier wise required skills
filepath_companies = "Companies"

new_dict_path = ''.join([json_folder,'/new_dict.json'])
file_dict_path = ''.join([json_folder,'/file_dict.json'])

if os.path.exists(new_dict_path):

    new_dict = read(new_dict_path)
    file_dict = read(file_dict_path)
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
                    if tmp in new_dict.keys():
                        new_dict[tmp].append(string)
                    else:
                        new_dict[tmp] = [string]

                    file_dict[string].append(tmp)

                c = c+1

    write(new_dict,new_dict_path)
    write(file_dict,file_dict_path)

# print file_dict

k = 0
carrier_dict={} # Stores Carrier that user can choose
skill_dict={}   # Stores Skills of user
with open('new_user.txt') as f:
    for line in f:
        if k==0:
            user_name = re.sub('\n','', line)
            if user_name not in skill_dict.keys():
                skill_dict[user_name] = []
        elif k==1:
            country = re.sub('\n','', line)
        elif k>=3:
            tmp = re.sub('\n','', line)
            if user_name in carrier_dict.keys():
                if tmp in new_dict.keys():
                    for i in range(len(new_dict[tmp])):
                        if new_dict[tmp][i] not in carrier_dict[user_name]:
                            carrier_dict[user_name].append(new_dict[tmp][i])
            else:
                if tmp in new_dict.keys():
                    carrier_dict[user_name] = new_dict[tmp]

            skill_dict[user_name].append(tmp)

        k = k+1

# print carrier_dict
# print skill_dict
username = user_name
final_dict = {}
car_final_dict = {}
os.system("clear")
num = []
for i in range(len(carrier_dict[username])):
    
    li = []
    cnt = 0
    for item in file_dict[carrier_dict[username][i]]:
        if item not in skill_dict[username]:
            li.append(item)
            cnt = cnt+1

    num.append(cnt)
    if cnt in final_dict.keys():
        final_dict[cnt].append(li)
        car_final_dict[cnt].append(carrier_dict[username][i])
    else:
        final_dict[cnt] = [li]
        car_final_dict[cnt] = [carrier_dict[username][i]]

num.sort()
cnt = 0
print "Hello",username,":-)"
print "\nBelow are most close carrier based on your current skills:"
for i in num:
    if cnt != 2:
        for j in range(len(car_final_dict[i])):
            print "\n-------------------------------------------------------------"
            print "If you want to become", car_final_dict[i][j], "you have to acquire following",str(i),"skills:\n"

            for k in range(len(final_dict[i][j])):
                print "    ->",final_dict[i][j][k]
            print "\n****Who will hire you?****"
            print "\nThe Companies that could hire you after becoming",car_final_dict[i][j],"in",country,"are:\n"
            dct = company_processing(filepath_companies,country,car_final_dict[i][j])
            for k in range(len(dct)):
                print "    ->",dct[k]
    else:
        quit()
    cnt = cnt+1
