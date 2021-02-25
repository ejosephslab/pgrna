#!/usr/bin/env python
# coding: utf-8

# Code for Cas13_pgRNA version 0.1

# In[2]:


# Reading the csv file for gRNA for sars cov2
import pandas as pd 
import re
import numpy as np
import os
import subprocess
# importing biopython
#from Bio import SeqIO


# In[3]:


import time
start = time.time()


# In[4]:


df = pd.read_csv("Input/CovidCasRxguides.csv")
df.head()


# In[5]:


# Checking the number of rows
len(df)


# In[6]:


## Finding the top quartile 
ng = len(df)# number of guide candidates

df = df.sort_values(by=["standardizedGuideScores"], ascending=False)


# In[7]:


n = ng//4 # To get the top quartile
df = df[:n] 


# In[8]:


# After sorting the length of dataframe
len(df)


# In[9]:


gRNA_list = list(df["GuideSeq"])
GuideName_list = list(df["MatchPos"])
score_list = list(df["standardizedGuideScores"])
# Function to compare two strings
def Str2MatchStr1(str1, str2):
    count = 0
    l = len(str1)
    for i in range(l):
        if str1[i]== str2[i]:
            count = count+1
    return count

match1 = []
match2 = []
similar_nt = []
score_gRNA_A = []
score_gRNA_B = []
GuideName_A = []
GuideName_B = []
l = len(gRNA_list)

for i in range(l-1):
    for j in range(i+1,l):
        #print(i,j)
            x = gRNA_list[i]
            row = gRNA_list[j]
            score1 = score_list[i]
            score2 = score_list[j]
            guideA = GuideName_list[i]
            guideB = GuideName_list[j]
            matchnumber = Str2MatchStr1(x, row) 
            if (matchnumber >= 17):
                match1.append(x)
                match2.append(row)
                similar_nt.append(matchnumber)
                score_gRNA_A.append(score1)
                score_gRNA_B.append(score2)
                GuideName_A.append(guideA)
                GuideName_B.append(guideB)
                #print('position', i,'gRNA', x, ' and ','position',j, 'match', row, 'similar nt', matchnumber)


# In[10]:


df_pairs = pd.DataFrame({'GuideA Position': GuideName_A,'Guide A':match1,'Guide A Score': score_gRNA_A,'GuideB Position': GuideName_B,'Guide B': match2, 'Guide B Score': score_gRNA_B,'Similar NT': similar_nt})
df_pairs.head()


# In[11]:


len(df_pairs)


# In[12]:


#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)


# In[13]:


## Doing Sequence Complement of Homologous Guide Sequences to find the Target position in SARS Cov2 genome
target_A_list = df_pairs["Guide A"].apply(revcom)
target_B_list = df_pairs["Guide B"].apply(revcom)

df_pairs["Target A"] = target_A_list
df_pairs["Target B"] = target_B_list


# In[14]:


df_pairs.head()


# In[15]:


def Mismatch_Alignment(a,b):
    mismatch = dif(a,b)
    for val in mismatch:
        list1 = list(a)
        list1[val] = 'X'
        a = ''.join(list1)

    return(a)


# In[16]:


def dif(a, b):
    return [i for i in range(len(a)) if a[i] != b[i]]


# In[17]:


Mismatch_position_list = []
for i in range(len(df_pairs)):
    a = df_pairs["Target A"][i]
    b = df_pairs["Target B"][i]
    Mismatch_position = Mismatch_Alignment(a,b)
    Mismatch_position_list.append(Mismatch_position)


# In[18]:


df_pairs["Mismatch_Position"] = Mismatch_position_list


# In[19]:


df_pairs.to_csv (r'df_pairs_sorted_top_quartile.csv', index = False, header=True)


# In[20]:


# Creating Decimal tobase 4 function in python
def is_zero(n):
    for d in n:
        if d != 0:
            return False
    return True

def modulo_div(n,original_base,destination_base):
    carry = 0
    for i in range(len(n)):
        d = n[i]
        d+=original_base*carry 
        carry = d%destination_base 
        d=(d//destination_base)
        n[i] = d
    return carry

def convertBase(n,original_base,destination_base):
    digits = []    
    while not is_zero(n):
        digits.insert(0,modulo_div(n,original_base,destination_base))
    return digits

def Dec2base4(inputNum):
    if inputNum <4:
        l1 = []
        num = inputNum//1
        l1.append(num)
    elif inputNum >= 4:
        n = str(inputNum)
        b = str(n)
        c = []

        for digit in b:
            c.append (int(digit))
        n = c
        original_base = 10
        destination_base = 4
        l1 = convertBase(n,original_base,destination_base)             
    return l1


# ### Design of optimised PgRNA

# #### GC Content
# Checking all the criteria for gRNA design
# GC content > 30% and <70%

# Defining a function
# GC content (or Guanine-Cytosine content)
def gc_content(genome):
    g_content = genome.count("G")
    c_content = genome.count("C")
    a_content = genome.count("A")
    t_content = genome.count("T")
    return 100*((g_content+c_content)/(a_content+c_content+g_content+t_content))


# #### CFD Scoring part

# In[22]:


df_mismatch = pd.read_csv("Input/Mismatch_Scores.csv")
#df_mismatch


# In[23]:


df_mismatch_scores = df_mismatch.set_index('Pair')
df_mismatch_scores


# In[24]:


df_adjacent_mismatch = pd.read_csv("Input/Adjacent_Mismatch.csv")
df_adjacent_mismatch


# In[26]:


## Finding the mismatch and calculating scores for single mismatch
def calc_cfd(wt,sg):
    M_1 = Find_single_Mismatch(wt,sg)
    #print("The position of Single Mismatch is",M_1)
    score = 1
    sg = revcom(sg)
    wt = revcom(wt)
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    #print("Reverse Complement sg:", sg)
    #print("sg:", sg)
    #print("DNA", wt)
    s_list = list(sg)
    wt_list = list(wt)
    #print(s_list)
    #print(wt_list)
    for i in M_1:
        j = len(sg)-(i+1)
        key = 't'+wt_list[j]+':s'+s_list[j]
        #print(key)
        position = str(i+1)
        #print(position)
        #print(df_mismatch_scores[position][key])
        score*= df_mismatch_scores[position][key]
    #print("Mismatch Penalty Score due to single Mismatch is:", score)
    return score


# In[27]:


def Find_Triple_Mismatch(a,b):
    a = a
    b = b
    ### Combining all
    def dif(a, b):
        return [i for i in range(len(a)) if a[i] != b[i]]
    mismatch = dif(a,b)
    for val in mismatch:
        list1 = list(a)
        list1[val] = 'X'
        a = ''.join(list1)
    import re
    l_3 = [m.start() for m in re.finditer('XXX', a)]
    return l_3


# In[28]:


def Find_Double_Mismatch(a,b):
    a = a
    b = b
    ### Combining all
    def dif(a, b):
        return [i for i in range(len(a)) if a[i] != b[i]]
    mismatch = dif(a,b)
    for val in mismatch:
        list1 = list(a)
        list1[val] = 'X'
        a = ''.join(list1)
    import re
    l_3 = [m.start() for m in re.finditer('XXX', a)]
    l_2 = [m.start() for m in re.finditer('XX', a)]
    def Diff(li1, li2): 
        li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2] 
        return li_dif
    l = Diff(l_2,l_3)
    return l


# In[29]:


def Find_single_Mismatch(a,b):
    a = a
    b = b
    ### Combining all
    def dif(a, b):
        return [i for i in range(len(a)) if a[i] != b[i]]
    mismatch = dif(a,b)
    for val in mismatch:
        list1 = list(a)
        list1[val] = 'X'
        a = ''.join(list1)
    import re
    l1_3 = [m.start() for m in re.finditer('XXX', a)]
    s1 = []
    for val in l1_3:
        #print(val)
        x = val+1
        y = val+2
        s = [x,y]
        s1.append(s)
    l1_3 = l1_3 + s1

    # function used for removing nested  
    # lists in python.
    output = [] 
    def reemovNestings(l):  
        for i in l: 
            if type(i) == list: 
                reemovNestings(i) 
            else: 
                output.append(i) 
        return output

    l1_3 = reemovNestings(l1_3) 

    l1_2 = [m.start() for m in re.finditer('XX', a)]
    s1 = []
    for val in l1_2:
        x = val+1
        s1.append(x)
    l1_2 = l1_2 + s1

    l1_consecutive = l1_3 + l1_2
    l1_consecutive = list(set(l1_consecutive))

    l1_1 = [m.start() for m in re.finditer('X', a)]


    def Diff(li1, li2): 
        li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2] 
        return li_dif


    l = Diff(l1_1,l1_consecutive)
    
    return l


# In[30]:


def Cal_Mismatch_Score(a,b):
    M_3 = Find_Triple_Mismatch(a,b)
    M_2 = Find_Double_Mismatch(a,b)
    M_1 = Find_single_Mismatch(a,b)

    score2 = 1
    score3 = 1
    score_1 = 1
    if (len(M_3) != 0) :
        for i in M_3:
            #print("The position of Triple mismatch is",i)
            j = 27-(i+1)
            score3 = score3*(df_adjacent_mismatch["Triple"][j])*(df_adjacent_mismatch["Triple"][j-1])*(df_adjacent_mismatch["Triple"][j-2])
    #print("The mismatch penalty due to Triple mismatch is",score3)
    if (len(M_2) != 0) :
        for i in M_2:
            #print("The position of Double mismatch is",i)
            k = 27-(i+1)
            score2 = score2*(df_adjacent_mismatch["Double"][k])*(df_adjacent_mismatch["Double"][k-1])
    #print("The mismatch penalty due to double mismatch is",score2)
    if (len(M_1) != 0) :
        #print(len(M_1))
        score_1 = calc_cfd(a,b)
        #print(score_1)
        #print("The mismatch penalty due to single mismatch is",score_1)
    #print(type(score_1))
    #print('score_1', score_1)
    Net_Score = float(score2)*float(score3)*float(score_1)
    return Net_Score


# ### Design of optimised PgRNA based on target activity at A and B greater than minimum of both columns Spacer A and B 

# In[31]:


## Design of optimised PgRNA Based on target activity at target A and B both > 0.75
l = len(df_pairs)
crRNA_list = []
Guide_score = []
Guide_scoreA = []
Guide_scoreB = []
for i in range(l):
    s1 = df_pairs["Guide A"][i]
    t1 = s1# storing as t1 because we are changing the s1 for the design of pgRNAs
    s2 = df_pairs["Guide B"][i]
    t2 = s2# storing as t2 because we are changing the s2 for the design of pgRNAs
    s1_score = df_pairs["Guide A Score"][i]
    s2_score = df_pairs["Guide B Score"][i]

    mismatch = dif(s1, s2)
    MMn = len(mismatch)
    #print(MMn)
    Com_MMn = ((4**MMn))
    BaseSwap_list = []
    for n in range(Com_MMn):
        #print(n)
        BaseSwap = Dec2base4(n)
        #print(BaseSwap);
        #print(type(BaseSwap))
        BaseSwap_list.append(BaseSwap)
        
    for val in BaseSwap_list:
        if(len(val) < MMn):
            zero_padding = MMn - len(val)
            listofzeros = [0]*zero_padding
            for n in listofzeros:
                val.insert(0, n)
            #print(val)
    #print("For Mismatch", MMn,"The list of replacer", BaseSwap_list)
    
    
    for pair in BaseSwap_list:
        l = len(pair)
        for j in range(l):
            if (pair[j] == 0):
                pair[j] = "A"
            elif (pair[j] == 1):
                pair[j] = "C"
            elif (pair[j] == 2):
                pair[j] = "G"
            elif (pair[j] == 3):
                pair[j] = "T"
        #print(pair)
        
    crRNA_list_1 = []
    for replacer in BaseSwap_list:
        i = 0
        for val in mismatch:
            #print(replacer[i])
            list1 = list(s1)
            list1[val] = replacer[i]
            i = i+1
            s1 = ''.join(list1)

        #print(s1)
        crRNA_list_1.append(s1)
        
        
    ## Checking the GC content of crRNA in CrRNA_list
    
    crRNA_list_1_corrected1 = []
    for i in range(len(crRNA_list_1)):
        gc_value = gc_content(crRNA_list_1[i])
        if (gc_value >30) and (gc_value < 70 ):
            crRNA_list_1_corrected1.append(crRNA_list_1[i])
            
            
    crRNA_list_1_corrected2 = []        
    # Criteria 2
    # To avoid repeats (≥4 consecutive ‘rU’ or ≥5 consecutive ‘rG’, ‘rC’, or ‘rA)
    l = len(crRNA_list_1_corrected1)
    for i in range(l):
        U_repeat = crRNA_list_1_corrected1[i].find('UUUU')
        G_repeat = crRNA_list_1_corrected1[i].find('GGGGG')
        C_repeat = crRNA_list_1_corrected1[i].find('CCCCC')
        A_repeat = crRNA_list_1_corrected1[i].find('AAAAA')
        if (U_repeat==-1) and (G_repeat==-1) and (C_repeat==-1) and (A_repeat==-1):
            crRNA_list_1_corrected2.append(crRNA_list_1_corrected1[i]) 
            
    #Checking the cfd cores
    Activity1 = min(df_pairs['Guide A Score'])  
    Activity2 = min(df_pairs['Guide B Score']) 
    crRNA_list_1_corrected3 = []
    Activity_guide = []
    Activity_guideA = []
    Activity_guideB = []
    for i in range(len(crRNA_list_1_corrected2)):
        mismatch_scoreA = Cal_Mismatch_Score(t1, crRNA_list_1_corrected2[i])
        Activity_scoreA  = s1_score* mismatch_scoreA
        mismatch_scoreB = Cal_Mismatch_Score(t2, crRNA_list_1_corrected2[i])
        Activity_scoreB  = s2_score* mismatch_scoreB
        
        Activity_guide_avg = (Activity_scoreA+Activity_scoreB)/2 # Taking the average of activity of Guide at target A and B
        
        if ( Activity_scoreA > Activity1) and (Activity_scoreB > Activity2):
            #print(mismatch_scoreA)
            #print(mismatch_scoreA)
            crRNA_list_1_corrected3.append(crRNA_list_1_corrected2[i])
            Activity_guide.append(Activity_guide_avg)
            Activity_guideA.append(Activity_scoreA)
            Activity_guideB.append(Activity_scoreB)
            
    ## Appending to final design RNA list
    d2 = {'Guide A' : t1, 'Guide B' : t2, 'Mismatch' : MMn , 'pgRNA List' : crRNA_list_1_corrected3}
    crRNA_list.append(d2)
    ## Appending the guide scores
    Guide_score.append(Activity_guide)
    Guide_scoreA.append(Activity_guideA)
    Guide_scoreB.append(Activity_guideB)
    
    


# In[32]:


df_crRNA = pd.DataFrame(crRNA_list)
df_crRNA


# In[33]:


x = 0
for i in range(len(df_crRNA)):
    #print(len(df_crRNA['pgRNA List'][i]))
    x = x + len(df_crRNA['pgRNA List'][i])
print(x)


# In[34]:


# Use list comprehension to convert the nested lists and get a full list of pgRNAs Activity Score
Guide_score_list = [ item for elem in Guide_score for item in elem]
len(Guide_score_list)


# In[35]:


# Use list comprehension to convert the nested lists and get a full list of pgRNAs Activity Score at Target A
Guide_scoreA_list = [ item for elem in Guide_scoreA for item in elem]
len(Guide_scoreA_list)


# In[36]:


# Use list comprehension to convert the nested lists and get a full list of pgRNAs Activity Score at Target A
Guide_scoreB_list = [ item for elem in Guide_scoreB for item in elem]
len(Guide_scoreB_list)


# In[37]:


df_crRNA = pd.DataFrame(crRNA_list)
df_crRNA


# In[38]:


#df_pairs["Guide A"] == df_crRNA["Guide A"]


# In[39]:


df_crRNA1 = df_crRNA
df_crRNA1['GuideA Position'] = df_pairs['GuideA Position']
df_crRNA1['Guide A Score'] = df_pairs['Guide A Score']
df_crRNA1['GuideB Position'] = df_pairs['GuideB Position']
df_crRNA1['Guide B Score'] = df_pairs['Guide B Score']
df_crRNA1 = df_crRNA1.rename(columns={"Guide A": "GuideA", "Guide B": "GuideB", "GuideA Position": "GuideAPosition", "GuideB Position": "GuideBPosition", "Guide A Score": "GuideAScore", "Guide B Score": "GuideBScore", "pgRNA List": "pgRNA"})
df_crRNA1


# In[40]:


df1 = df_crRNA1[['GuideA', 'pgRNA']].copy()
GuideAs = []
pgRNAs = []
for _, row in df1.iterrows():
    GuideA = row.GuideA
    
    for val in row.pgRNA:
        GuideAs.append(GuideA)
        pgRNAs.append(val)
df2 = pd.DataFrame({"GuideA": GuideAs,
    "pgRNA": pgRNAs
})
df2.head(10)


# In[41]:


df1 = df_crRNA1[['GuideB', 'pgRNA']].copy()
GuideBs = []
pgRNAs = []
for _, row in df1.iterrows():
    GuideB = row.GuideB
    
    for val in row.pgRNA:
        GuideBs.append(GuideB)
        pgRNAs.append(val)
df2 = pd.DataFrame({"GuideB": GuideBs,
    "pgRNA": pgRNAs
})
df2.head(10)


# In[42]:


df1 = df_crRNA1[['GuideAPosition', 'pgRNA']].copy()
GuideAPositions = []
pgRNAs = []
for _, row in df1.iterrows():
    GuideAPosition = row.GuideAPosition
    
    for val in row.pgRNA:
        GuideAPositions.append(GuideAPosition)
        pgRNAs.append(val)
df2 = pd.DataFrame({"GuideAPosition": GuideAPositions,
    "pgRNA": pgRNAs
})
df2.head(10)


# In[43]:


df1 = df_crRNA1[['GuideBPosition', 'pgRNA']].copy()
GuideBPositions = []
pgRNAs = []
for _, row in df1.iterrows():
    GuideBPosition = row.GuideBPosition
    
    for val in row.pgRNA:
        GuideBPositions.append(GuideBPosition)
        pgRNAs.append(val)
df2 = pd.DataFrame({"GuideBPosition": GuideBPositions,
    "pgRNA": pgRNAs
})
df2.head(10)


# In[44]:


df1 = df_crRNA1[['GuideAScore', 'pgRNA']].copy()
GuideAScores = []
pgRNAs = []
for _, row in df1.iterrows():
    GuideAScore = row.GuideAScore
    
    for val in row.pgRNA:
        GuideAScores.append(GuideAScore)
        pgRNAs.append(val)
df2 = pd.DataFrame({"GuideAScore": GuideAScores,
    "pgRNA": pgRNAs
})
df2.head(10)


# In[45]:


df1 = df_crRNA1[['GuideBScore', 'pgRNA']].copy()
GuideBScores = []
pgRNAs = []
for _, row in df1.iterrows():
    GuideBScore = row.GuideBScore
    
    for val in row.pgRNA:
        GuideBScores.append(GuideBScore)
        pgRNAs.append(val)
df_crRNA_expanded = pd.DataFrame({"GuideA": GuideAs, "GuideB": GuideBs,
    "pgRNA": pgRNAs
})
df_crRNA_expanded.head(10)


# In[46]:


df_crRNA_expanded["GuideScore"] = Guide_score_list
df_crRNA_expanded["GuideScore_At_TargetA"] = Guide_scoreA_list
df_crRNA_expanded["GuideScore_At_TargetB"] = Guide_scoreB_list
df_crRNA_expanded["GuideAScore"] = GuideAScores
df_crRNA_expanded["GuideBScore"] = GuideBScores
df_crRNA_expanded["GuideAPosition"] = GuideAPositions
df_crRNA_expanded["GuideBPosition"] = GuideBPositions


# In[47]:


# deleting df_crRNA1 dataframe
del df_crRNA1


# In[48]:


df_crRNA1 = df_crRNA_expanded


# In[49]:


### Adding a column for the GC content
GC_content_list = df_crRNA1["pgRNA"].apply(gc_content)


# In[50]:


df_crRNA1["GC_content%"] = GC_content_list


# In[51]:


df_crRNA1


# #### Doing Multiple sequences call in RNAfold

# In[52]:


# Initializing list  
list1 = df_crRNA1['pgRNA'].values.tolist()
  
len(list1)


# In[53]:


# Criteria 3
# Min Free energy of hainpin + pgRNA > -5
##Adding the hairpin structure to the pgRNAs in the final list2
list3 = []
for val in list1:
    Hainpn = "AACCCCTACCAACTGGTCGGGGTTTGAAAC"
    crRNA_combined = Hainpn + val
    list3.append(crRNA_combined)


# In[54]:


# Converting the list3 to a text file
with open('crRNA_combined.txt', 'w') as f:
    for item in list3:
        f.write("%s\n" % item)


# In[55]:


# Converting the text file to a fasta file

fileInput = open("crRNA_combined.txt", "r")
fileOutput = open("out_crRNA_combined.fsa", "w")
#Seq count
count = 1 

#Loop through each line in the input file
print("Converting to FASTA...")
for strLine in fileInput:

    #Strip the endline character from each input line
    strLine = strLine.rstrip("\n")

    #Output the header
    fileOutput.write(">" + str(count) + "\n")
    fileOutput.write(strLine + "\n")

    count = count + 1
print ("Done.")

#Close the input and output file
fileInput.close()
fileOutput.close()


# In[56]:


### Function to calculate the Minimum Free energy 
def cal_Min_Eng(x):
    s = x.split("\n")[0]
    s_negative = s.find("-")
    numeric_filter = filter(str.isdigit, s)
    numeric_string = "".join(numeric_filter)
    l = list(numeric_string)
    length_l = len(l)
    l.insert(length_l-2, ".")
    Mng = ''.join(l)
    Mng = float(Mng)
    if (s_negative == -1):
        Mng = Mng
    else:
        Mng = -Mng
        
    return Mng


# In[57]:


### Function to calculate the Secondary Structure Using RNAFold Software
def cal_Structure(s):
    s1 = s.split(" ")
    Structure = s1[0]
    return Structure


# In[58]:


command1 = "RNAfold < out_crRNA_combined.fsa >> energyoutput.txt"
subprocess.call(command1, shell=True) 


# In[59]:


command2 = "del *_ss.ps"
subprocess.call(command2, shell=True) 


# In[60]:


f=open("energyoutput.txt", "r")
contents =f.readlines()
l = len(contents)
sequence_count = []
energy_count = []
for i in range(l):
    For_Query = contents[i].find(">")
    #print(For_Query)
    if (For_Query == 0):
        sequence_count.append(contents[i+1])
        energy_count.append(contents[i+2])
        
f.close()


# In[61]:


energy_list = []
Structure_list = []
for val in energy_count:
    energy = cal_Min_Eng(val)
    energy_list.append(energy)
    Structure = cal_Structure(val)
    Structure_list.append(Structure)


# In[62]:


df_crRNA1["Free Energy"] = energy_list
df_crRNA1["Secondary Structure"] = Structure_list


# In[63]:


HairpinNRG = -12.8
df_crRNA2 = df_crRNA1.loc[(df_crRNA1['Free Energy']- HairpinNRG) > -5]
len( df_crRNA2)


# In[64]:


df_crRNA2 = df_crRNA2.reset_index(drop=True)
df_crRNA2


# In[65]:


df_crRNA2.to_csv (r'df_crRNA_withEnergy_filtered.csv', index = False, header=True)


# In[66]:


list_filtered = df_crRNA2['pgRNA'].values.tolist()
  
len(list_filtered)


# In[67]:


reverse_complement_rna_list = []
for val in list_filtered:
    # Taking the reverse complement of pg-DNA sequence
    s = revcom(val)
    s = s.replace('T', 'U')
    reverse_complement_rna_list.append(s)
    
    


# In[68]:


len(reverse_complement_rna_list)


# In[69]:


# Converting the reverse_complement_rna_list into a text file
with open('pgRNA_filtered.txt', 'w') as f:
    for item in reverse_complement_rna_list:
        f.write("%s\n" % item)


# In[70]:


# Converting the text file to a fasta file

fileInput = open("pgRNA_filtered.txt", "r")
fileOutput = open("out_pgRNA_filtered.fsa", "w")
#Seq count
count = 1 

#Loop through each line in the input file
print("Converting to FASTA...")
for strLine in fileInput:

    #Strip the endline character from each input line
    strLine = strLine.rstrip("\n")

    #Output the header
    fileOutput.write(">" + str(count) + "\n")
    fileOutput.write(strLine + "\n")

    count = count + 1
print ("Done.")

#Close the input and output file
fileInput.close()
fileOutput.close()


# ### Blast search

# In[71]:


command = 'blastn -query out_pgRNA_filtered.fsa -task blastn-short -db C:/Users/18064/Documents/Blast/db/Human_NCBI_rnadb -outfmt "7 qacc sacc qstart qend sstart send" -out pgRNA_NCBI_blast.out'
subprocess.call(command, shell=True)


# In[72]:


f=open('pgRNA_NCBI_blast.out', "r")
contents =f.readlines()
query_list =[]
hits_list = []
for x in contents:
    For_Query = x.find("Query")
    if (For_Query != -1):
        query_list.append(x)
        #print("For Sequence",x)
    hits_count = x.find("hits found")
    if (hits_count != -1):
        hits_list.append(x)
        #print(x)

# Making a list of query sequence
new_query_list =[]
for val in query_list:
    x = val.split("\n")[0].split("#")[1].split(":")[1].split(" ")[1]
    x = int(x)
    new_query_list.append(x)

# Making a list of blast hits in host genome for the query sequences
new_hit_list =[]
for val in hits_list:
    x = val.split("\n")[0].split(" ")[1]
    x = int(x)
    new_hit_list.append(x)


# In[73]:


blast_hits = pd.DataFrame({'GuideNumber': new_query_list,'Hits': new_hit_list})
blast_hits.head()


# In[74]:


len(blast_hits)


# In[75]:


df_crRNA2["HumanBlastHits"] = new_hit_list


# In[76]:


df_crRNA2.to_csv (r'df_crRNA_withhits.csv', index = False, header=True)


# In[77]:


f=open("pgRNA_NCBI_blast.out", "r")
contents =f.readlines()
Blast_list = []
for i in range(len(contents)):
    target = contents[i].find("#")
    if (target == -1):
        Blast_list.append(i)


# In[78]:


seq_list = []
accession_list = []
qstart_list = []
qend_list = []
sstart_list = []
send_list = []
for val in Blast_list:
    line = contents[val]
    seq = line.split('\t')[0]
    accession = line.split('\t')[1]
    qstart = line.split('\t')[2]
    qend = line.split('\t')[3]
    sstart = line.split('\t')[4]
    send = line.split('\t')[5].split('\n')[0]
    seq_list.append(seq)
    accession_list.append(accession)
    qstart_list.append(qstart)
    qend_list.append(qend)
    sstart_list.append(sstart)
    send_list.append(send)
    
d = {'Sequence':seq_list,'Accession':accession_list,'qstart':qstart_list,'qend':qend_list,'sstart':sstart_list,'send':send_list}
df_blast = pd.DataFrame(d)
df_blast.head()


# In[79]:


# This is for Cas13d
range_list = []
seq_4_range_list = []
l = len(sstart_list)
for i in range(l):
    if (sstart_list[i] < send_list[i]): # For RNA, hit same strand
        #print('strand = plus')

        range_blast = accession_list[i] + '  '+ str(int(sstart_list[i])-int(qstart_list[i])+1) + '-' + str(int(send_list[i])+ 23 - int(qend_list[i])) + '  '+ 'plus'
        range_list.append(range_blast)
        seq_4_range_list.append(seq_list[i])


# In[80]:


range_list


# In[81]:


len(range_list)


# In[82]:


#seq_4_range_list


# In[83]:


# Converting the list to a text file
with open('range_NCBIBlasthits.txt', 'w') as f:
    for item in range_list:
        f.write("%s\n" % item)


# In[84]:


command1 = 'blastdbcmd -db C:/Users/18064/Documents/Blast/db/Human_NCBI_rnadb -entry_batch range_NCBIBlasthits.txt -outfmt "%s" -out NCBI_Human_RNA_Range.txt'
subprocess.call(command1, shell=True)


# In[85]:


## Reading the NCBI_Human_RNA_Range.txt file
f_rangeRNA=open("NCBI_Human_RNA_Range.txt", "r")
contents_rangeRNA =f_rangeRNA.readlines()
Human_rangeRNA_list = []
for val in contents_rangeRNA:
    val1 = (val.split('\n')[0])[0:23]
    Human_rangeRNA_list.append(val1)   


# In[86]:


len(Human_rangeRNA_list)


# In[87]:


seq_blast_cfd = []
GuideScore_blast = []
for val in seq_4_range_list:
    i = (int(val) - 1)
    #print(i)
    seq_blast = revcom(df_crRNA2.iloc[i].loc['pgRNA'])
    guide_score = df_crRNA2.iloc[i].loc['GuideScore']
    seq_blast_cfd.append(seq_blast)
    GuideScore_blast.append(guide_score)


# In[88]:


Blast_cfd_df = pd.DataFrame({'GuideNumber': seq_4_range_list,'SequenceComplementGuide': seq_blast_cfd,'HumanRNA_range': Human_rangeRNA_list,'GuideActivity': GuideScore_blast})


# In[89]:


Blast_cfd_df.head()
    


# In[91]:


## Calculating the cfd score for the host transcriptome region corresponding to blast hits
x = len(Human_rangeRNA_list)
cfd_humanBlast_list = []
Activity_humanBlast_list = []
for i in range(x):
    #print(i)
    try:
        cfd_humanBlast = Cal_Mismatch_Score(Blast_cfd_df["SequenceComplementGuide"][i], Blast_cfd_df["HumanRNA_range"][i])
        Activity_humansite = cfd_humanBlast*Blast_cfd_df["GuideActivity"][i]
    except IndexError:
        cfd_humanBlast = 'null'
        Activity_humansite = "null"
    cfd_humanBlast_list.append(cfd_humanBlast)
    Activity_humanBlast_list.append(Activity_humansite)


# In[92]:


Blast_cfd_df["CFD_Score"] = cfd_humanBlast_list
Blast_cfd_df["ActivityScore"] = Activity_humanBlast_list


# In[93]:


Blast_cfd_df.head()


# In[94]:


Blast_cfd_df["Guide"] = Blast_cfd_df["SequenceComplementGuide"].apply(revcom)


# In[95]:


Blast_cfd_df.head()


# In[96]:


Blast_cfd_df.to_csv (r'BlastHumanRnaHits_Activity.csv', index = False, header=True)


# In[97]:


end = time.time()
print(end - start)


# In[ ]:





# In[ ]:





# In[ ]:




