
import math
import os
import zlib
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from scipy.optimize import fsolve
import copy

class Node:
    def __init__(self, data,counting):
        self.str1=data
        self.counting=counting

#for k = 10 in SI
def calcExactDistFromSi(t, synIn):
  if (synIn == 0.0):
               print ("fgfgg")
               return(1000.0)
  #if (synIn == 1):
  #        return 0
  else:
    return(math.exp(2*t)*synIn - (10*t**18 + 90*t**17 + 645*t**16 + 2580*t**15 + 8604*t**14 + 20076*t**13 + 39804*t**12 +
  59706*t**11 + 76502*t**10 + 76502*t**9 + 64998*t**8 + 43332*t**7 + 24108*t**6 + 10332*t**5 + 3552*t**4 + 888*t**3 +
  162*t**2 + 18*t + 1)/(t + 1)**19)


def SyntenyIndex_indels(Genom1, Genom2, k):
    total = 0
    t = 0
    for item in range(k, len(Genom1) - k):
        if (Genom1[item] in Genom2):
          total += 1
          set1 = set(Genom1[item - k:item])
          set2 = set(Genom1[item + 1:item + k + 1])
          set3 = set1.union(set2)
          item1 = Genom2.index(Genom1[item])
          set12 = set(Genom2[item1 - k:item1])
          set22 = set(Genom2[item1 + 1:item1 + k + 1])
          set32 = set12.union(set22)
          Intersection = set3.intersection(set32)
          lenintersection = len(Intersection)
          b = (1 / (2 * k))
          JSI = b * lenintersection
          t += JSI
    if (total == 0):
      return 0,0
    syntenyindex = (1 / total) * t
    return syntenyindex, total

   
def calc_lz_distance(source,seq,flagGen,start):
    count=start
    if flagGen==0:
        i=1
        count=count+1
    else:
        seq=source+seq
        i=len(source)+1
   
    subStr=seq[i]
    while i<len(seq):
        j=i
        subStr=seq[j]
        while subStr in seq[0:j]:
            i=i+1
            if i>len(seq):
                break
            subStr=seq[j:i+1]
        count=count+1
        i=i+1
               
    return count
    #i=0
    #while i<len(seq):
    #    subStr=seq[i]
    #    if subStr in dic and dic[subStr]==0:
    #        dic[subStr]=1
    #        i=i+1
    #        continue
    #    j=i
    #    if flagGen==0:
    #        while subStr in seq[0:j] or subStr in dic.keys():
    #            i=i+1
    #            if i>len(seq):
    #                break
    #            subStr=seq[j:i+1]
    #        dic[subStr]=1
    #        i=i+1
    #    else:
    #        while subStr in source or subStr in dic.keys():
    #            i=i+1
    #            if i>len(seq):
    #                break
    #            subStr=seq[j:i+1]
    #        dic[subStr]=1
    #        i=i+1
    #count=sum(value==1 for value in dic.values())
   

def lz_distance(gen1,gen2,s,q):
    sq=calc_lz_distance(gen1,gen2,1,s)
    qs=calc_lz_distance(gen2,gen1,1,q)
   
   
    #d=max(sq-s,qs-q)/max(sq,qs)
   
    #d=max(sq-s,qs-q)/max(s,q)
    d=min(((sq-s)/q),((qs-q)/s))
    return d




def calc_ncd(genom1,genom2):
    compressed_str1=zlib.compress(genom1.encode())
    compressed_str2=zlib.compress(genom2.encode())

    len1=len(compressed_str1)
    len2=len(compressed_str2)

    combined_len1=len(zlib.compress(genom1.encode()+genom2.encode()))
    combined_len2=len(zlib.compress(genom2.encode()+genom1.encode()))
    ncd=(combined_len1+combined_len2-2*min(len1,len2))/(2*max(len1,len2))
    return ncd

def calculate_len_for_SI (arr, k,names):
    flag=0
   
    siLength = [[0]*len(arr) for i in range(len(arr))]  
    fileForNames=open("names.txt","w")
    for name in range(len(names)):
       fileForNames.write(names[name]+"\n")
    fileForNames.close()
    for x in range (len(arr)):
        for y in range (x, len(arr)):
            if (x == y):
                siLength[x][x] = 0
            else:
                ncd=np.zeros((len(arr[x]), len(arr[y])))
                for i in range(len(arr[x])):
                    print("i=",i," x=",x," y=",y)
                    for j in range(len(arr[y])):
                        ncd[i][j]=calc_ncd(arr[x][i],arr[y][j])
                       # lz[i][j]=lz_distance(arr[x][i].str1,arr[y][j].str1,arr[x][i].counting,arr[y][j].counting)
                ncd1 = csr_matrix(ncd)
                minWeight = min_weight_full_bipartite_matching(ncd1)[1]
                if(len(arr[x])<len(arr[y])):
                    genomOrder1=np.arange(len(arr[x]))
                    genomOrder2=np.ones(len(arr[y]))*(-1)
                else:
                    genomOrder1=np.arange(len(arr[y]))
                    genomOrder2=np.ones(len(arr[x]))*(-1)
                for i in range(len(minWeight)):
                    #if(len(arr[x][i])<len(arr[y][minWeight[i]])):
                    #    smaller=arr[x][i]
                    #    bigger=arr[y][minWeight[i]]
                    #else:
                    #    smaller=arr[y][minWeight[i]]
                    #    bigger=arr[x][i]
                    #cutOff=0.1*(len(smaller)/len(bigger))*1.2+0.6*((len(bigger)-len(smaller))/len(bigger))
                    #if(ncd[i][minWeight[i]]<cutOff):
                    genomOrder2[minWeight[i]]=i
                result, total = SyntenyIndex_indels(list(genomOrder1), list(genomOrder2), k)
                if (result == 0.0):
                    dist = []
                    dist.append(1000)
                elif (result == 1):
                    dist = []
                    dist.append(0)
                else:
                    if (k==10):  
                        dist = fsolve(calcExactDistFromSi, 1.0, args = result)  #1.0 is the initial value for t
#             elif (k == 6):
#                dist = fsolve(calcExactDistFromSi_6, 1.0, args = result)  #1.0 is the initial value for t
#             else:
#               dist = fsolve(calcforsi, 1.0, args = result) #1.0 is the initial value for t  
                siLength[x][y] = round(dist[0],6)
                siLength[y][x] = round(dist[0],6)    
                if(flag==0):
                    fileForDistances=open("Distances2.txt","w")
                    fileForDistances.write(str(siLength[x][y])+ "\t")
                    flag=1
                    fileForDistances.close()
                else:
                    fileForDistances=open("Distances2.txt","a")
                    fileForDistances.write(str(siLength[x][y])+ "\t")
                    fileForDistances.close()
        fileForDistances=open("Distances2.txt","a")
        fileForDistances.write("\n")
        fileForDistances.close()
    newStr = str(len(arr)) + "\n"      
    for x in range (len(names)):
       newStr += names[x]
       for y in range (len(arr)):
            newStr+= " " + str(siLength[x][y])  
       newStr+="\n"    
    f = open("infile", "w")
    f.write(newStr)
    f.close()  
    return siLength, newStr

def createFilesForTreedist (n, scaleexp ,sizeofGenom,k,leaves,names):
   
   #call for function that creats matrix in "infile"
   calculate_len_for_SI(leaves,k,names)
  #call neighbor  
   os.system("echo y| phylip neighbor")
  #rename files
   old_name = r"outfile"
   new_name = r"new_tree" +".txt"
   os.rename(old_name, new_name)
   old_name = r"outtree"
   new_name =  r"intree2"
   os.rename(old_name, new_name)
  #call treedist
   stat = os.system("phylip treedist < phylip-param.txt > dist.out.txt");
  #rename files
   old_name = r"outfile"
   new_name = r"treedist" + ".txt"
   os.rename(old_name, new_name)
   old_name = r"intree"
   new_name = r"first"+".txt"
   os.rename(old_name, new_name)
   old_name = r"intree2"
   new_name = r"second"+".txt"
   os.rename(old_name, new_name)
   return

#def calculate_distances (n, scaleexp ,sizeofGenom,k,root,leaves):
#    #creates tree, calls to neighbor and treedist
#    tree = createFilesForTreedist(n, scaleexp ,sizeofGenom,k,root,leaves)
#    #open the result of treedist
#    temp = open("treedist" +".txt",'r')
#    f = temp.read()
    #get rf
#    RF = f.lstrip('1')
#    print (RF)
#    first = open ("first"+ ".txt", 'r')
#    first = first.read()
#    second = open ("second" + ".txt", 'r')
#    second = second.read()
    #count1 = first.count("(") - 1 - binaryTree(first)
#    count1 = first.count("(")
#    count2 = second.count ("(")  
#    common_edges = (count1 + count2 - int(RF))/2
#    result = int(RF)/ (count1 + count2)
#    return result,RF


n=30
folder_path="/home/user/Desktop/sapir/genoms/"

translate_table=str.maketrans('acgt','tgca')
matrix_genes=[]
matrix_names=[]
file_list=os.listdir(folder_path)
for file_name in file_list:
    genes=[]
    names=[]
    file_path=os.path.join(folder_path,file_name)
    file1=open(file_path,'r')
    content=file1.read()
    content=content.split('\n')
    genes=content[0].split(',')
    names=content[1].split(',')
    for i in range(len(names)):
        if('(C)' in names[i]):
            temp=genes[i]
            print(len(temp))
            temp=temp.translate(translate_table)
            temp=temp[::-1]
            #s=calc_lz_distance("", temp, 0,0)
            genes[i]=temp
        else:
            s=calc_lz_distance("", genes[i], 0,0)
            genes[i]=genes[i]
    file1.close()        
    matrix_genes.append(genes)
    matrix_names.append(names)
    names=[]
for file_name in file_list:
    name="species_"+file_name[0:-4]
    names.append(name)  
tree = createFilesForTreedist(n, 0.1 ,100,10,matrix_genes,names)