#!/usr/bin/env python
# coding: utf-8

# In[14]:


# Reading the csv file for gRNA for HPV16
import pandas as pd 
import re
import os
import subprocess
import pickle


# In[15]:


from Bio import SeqIO
for sequence in SeqIO.parse('HPV16.fasta', "fasta"):
    #print(sequence.seq)
    print(len(sequence),'nuceoliotides')


# In[16]:


df1 = pd.read_csv("HIV16_scores.csv")
df1.head()


# In[17]:


len(df1)


# In[18]:


df = df1[['sgRNA Cut Position (1-based)', 'Orientation', 'sgRNA Sequence', 'On-Target Efficacy Score']].copy()


# In[19]:


df.head()


# In[20]:


## Finding the top quartile 
ng = len(df)# number of guide candidates

df = df.sort_values(by=["On-Target Efficacy Score"], ascending=False)


# In[21]:


gRNA_list = list(df["sgRNA Sequence"])
GuideName_list = list(df["sgRNA Cut Position (1-based)"])
score_list = list(df["On-Target Efficacy Score"])
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
            if (matchnumber >= 15):
                match1.append(x)
                match2.append(row)
                similar_nt.append(matchnumber)
                score_gRNA_A.append(score1)
                score_gRNA_B.append(score2)
                GuideName_A.append(guideA)
                GuideName_B.append(guideB)
                #print('position', i,'gRNA', x, ' and ','position',j, 'match', row, 'similar nt', matchnumber)


# In[22]:


df_pairs = pd.DataFrame({'GuideA Position': GuideName_A,'Guide A':match1,'Guide A Score': score_gRNA_A,'GuideB Position': GuideName_B,'Guide B': match2, 'Guide B Score': score_gRNA_B,'Similar NT': similar_nt})
df_pairs


# In[23]:


len(df_pairs)


# In[24]:


#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)


# In[25]:


## Doing Sequence Complement of Homologous Guide Sequences to find the Target position in SARS Cov2 genome
target_A_list = df_pairs["Guide A"].apply(revcom)
target_B_list = df_pairs["Guide B"].apply(revcom)

df_pairs["Target A"] = target_A_list
df_pairs["Target B"] = target_B_list


# In[26]:


df_pairs.head()


# In[27]:


def Mismatch_Alignment(a,b):
    mismatch = dif(a,b)
    for val in mismatch:
        list1 = list(a)
        list1[val] = 'X'
        a = ''.join(list1)

    return(a)


# In[28]:


def dif(a, b):
    return [i for i in range(len(a)) if a[i] != b[i]]


# In[29]:


Mismatch_position_list = []
for i in range(len(df_pairs)):
    a = df_pairs["Target A"][i]
    b = df_pairs["Target B"][i]
    Mismatch_position = Mismatch_Alignment(a,b)
    Mismatch_position_list.append(Mismatch_position)


# In[30]:


df_pairs["Mismatch_Position"] = Mismatch_position_list


# In[31]:


df_pairs.to_csv (r'df_pairs_sorted_top_quartile.csv', index = False, header=True)


# In[32]:


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

# In[33]:


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

# In[34]:


# df_mismatch = pd.read_csv("Mismatch_Scores.csv")
# df_mismatch_scores = df_mismatch.set_index('Pair')
# df_mismatch_scores


# In[35]:


with open('mismatch_score.pkl', "rb") as f:
    lines = [line.rstrip(b"\r\n") for line in f.readlines()]
    mm_scores = pickle.loads(b"\n".join(lines))


# #### Reading the PAM scores from the pickle file

# In[36]:


with open('pam_scores.pkl', "rb") as f:
    lines = [line.rstrip(b"\r\n") for line in f.readlines()]
    pam_scores = pickle.loads(b"\n".join(lines))


# In[37]:


print(pam_scores)


# In[39]:


#Calculates CFD score
def calc_cfd(wt,sg,pam):
    #mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    #print ("sg:", sg)
    #print ("DNA", wt)
    s_list = list(sg)
    #print("s_list",s_list)
    wt_list = list(wt)
    #print("wt_list",wt_list)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            #print(key)
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score)


# #### Final Function to calculate the cfd score of Cas9 given two input sequence, polyvalent gRNA (20mer +PAM) and Spacer gDNA (20mer +PAM)

# In[40]:


# a is the Spacer gDNA/RNA and b is the pgDNA
def Cal_Mismatch_Score(a,b):
    #wt = a[:-3] 
    wt = a
    #off = b
    #sg = off[:-3]
    sg = b
    #pam = off[-2:] # Assuming pam as NGG
    pam = 'GG'
    score = calc_cfd(wt,sg,pam)
    return (score)


# In[41]:


Cal_Mismatch_Score('GGGAAAAAAAAAAAAAAAAAA','AAAAAAAAAAAAAAAAAAAAA')


# In[42]:


Cal_Mismatch_Score('CAGTACAGTAGTGGAAGTGG','CAGTACAGTGGTGGAAGTGG')


# In[43]:


Cal_Mismatch_Score('CAGTACAGTAGTGGAAGTGG','CAGTACAGTAGTGGAAGTGG')


# In[44]:


len('ATTAGCAGAACTACACACCA')


# #### Design of optimised PgRNA based on target activity at A and B greater than minimum of both columns Spacer A and B 

# In[45]:


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
        #print(gc_value)
        if (gc_value >30) and (gc_value < 70 ):
            #print(gc_value)
            crRNA_list_1_corrected1.append(crRNA_list_1[i])
            #print(crRNA_list_1_corrected1)
            
            
    crRNA_list_1_corrected2 = []        
    # Criteria 2
    # To avoid repeats (≥4 consecutive ‘rU’ or ≥5 consecutive ‘rG’, ‘rC’, or ‘rA)
    l1 = len(crRNA_list_1_corrected1)
    for i in range(l1):
        U_repeat = crRNA_list_1_corrected1[i].find('UUUU')
        G_repeat = crRNA_list_1_corrected1[i].find('GGGGG')
        C_repeat = crRNA_list_1_corrected1[i].find('CCCCC')
        A_repeat = crRNA_list_1_corrected1[i].find('AAAAA')
        if (U_repeat==-1) and (G_repeat==-1) and (C_repeat==-1) and (A_repeat==-1):
            crRNA_list_1_corrected2.append(crRNA_list_1_corrected1[i])
    #print(crRNA_list_1_corrected2)
            
    #Checking the cfd cores
    Activity1 = min(df_pairs['Guide A Score'])  
    Activity2 = min(df_pairs['Guide B Score']) 
    crRNA_list_1_corrected3 = []
    Activity_guide = []
    Activity_guideA = []
    Activity_guideB = []
    for i in range(len(crRNA_list_1_corrected2)):
        mismatch_scoreA = Cal_Mismatch_Score(t1, crRNA_list_1_corrected2[i])
        #print(mismatch_scoreA)
        Activity_scoreA  = s1_score* mismatch_scoreA
        mismatch_scoreB = Cal_Mismatch_Score(t2, crRNA_list_1_corrected2[i])
        #print(mismatch_scoreA)
        Activity_scoreB  = s2_score* mismatch_scoreB
        
        Activity_guide_avg = (Activity_scoreA+Activity_scoreB)/2 # Taking the average of activity of Guide at target A and B
        
        if ( Activity_scoreA > Activity1) and (Activity_scoreB > Activity2):
            #print(mismatch_scoreA)
            #print(mismatch_scoreA)
            crRNA_list_1_corrected3.append(crRNA_list_1_corrected2[i])
            Activity_guide.append(Activity_guide_avg)
            Activity_guideA.append(Activity_scoreA)
            Activity_guideB.append(Activity_scoreB)
    #print(crRNA_list_1_corrected3)
            
    ## Appending to final design RNA list
    d2 = {'Guide A' : t1, 'Guide B' : t2, 'Mismatch' : MMn , 'pgRNA List' :crRNA_list_1_corrected3 }
    crRNA_list.append(d2)
    ## Appending the guide scores
    Guide_score.append(Activity_guide)
    Guide_scoreA.append(Activity_guideA)
    Guide_scoreB.append(Activity_guideB)
    
    


# In[46]:


df_crRNA = pd.DataFrame(crRNA_list)
df_crRNA


# In[47]:


x = 0
for i in range(len(df_crRNA)):
    #print(len(df_crRNA['pgRNA List'][i]))
    x = x + len(df_crRNA['pgRNA List'][i])
print(x)


# In[48]:


# Use list comprehension to convert the nested lists and get a full list of pgRNAs Activity Score
Guide_score_list = [ item for elem in Guide_score for item in elem]
len(Guide_score_list)


# In[49]:


# Use list comprehension to convert the nested lists and get a full list of pgRNAs Activity Score at Target A
Guide_scoreA_list = [ item for elem in Guide_scoreA for item in elem]
len(Guide_scoreA_list)


# In[50]:


# Use list comprehension to convert the nested lists and get a full list of pgRNAs Activity Score at Target A
Guide_scoreB_list = [ item for elem in Guide_scoreB for item in elem]
len(Guide_scoreB_list)


# In[51]:


#df_pairs["Guide A"] == df_crRNA["Guide A"]


# In[52]:


df_crRNA1 = df_crRNA
df_crRNA1['GuideA Position'] = df_pairs['GuideA Position']
df_crRNA1['Guide A Score'] = df_pairs['Guide A Score']
df_crRNA1['GuideB Position'] = df_pairs['GuideB Position']
df_crRNA1['Guide B Score'] = df_pairs['Guide B Score']
df_crRNA1 = df_crRNA1.rename(columns={"Guide A": "GuideA", "Guide B": "GuideB", "GuideA Position": "GuideAPosition", "GuideB Position": "GuideBPosition", "Guide A Score": "GuideAScore", "Guide B Score": "GuideBScore", "pgRNA List": "pgRNA"})
df_crRNA1


# In[53]:


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
df2


# In[54]:


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
df2


# In[55]:


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
df2


# In[56]:


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
df2


# In[57]:


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
df2


# In[58]:


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
df_crRNA_expanded


# In[59]:


df_crRNA_expanded["GuideScore"] = Guide_score_list
df_crRNA_expanded["GuideScore_At_TargetA"] = Guide_scoreA_list
df_crRNA_expanded["GuideScore_At_TargetB"] = Guide_scoreB_list
df_crRNA_expanded["GuideAScore"] = GuideAScores
df_crRNA_expanded["GuideBScore"] = GuideBScores
df_crRNA_expanded["GuideAPosition"] = GuideAPositions
df_crRNA_expanded["GuideBPosition"] = GuideBPositions


# In[60]:


# deleting df_crRNA1 dataframe
del df_crRNA1


# In[61]:


df_crRNA1 = df_crRNA_expanded


# In[62]:


### Adding a column for the GC content
GC_content_list = df_crRNA1["pgRNA"].apply(gc_content)
df_crRNA1["GC_content%"] = GC_content_list
df_crRNA1


# #### Doing Multiple sequences call in RNAfold

# In[63]:


# Initializing list  
list1 = df_crRNA1['pgRNA'].values.tolist()
  
len(list1)


# In[64]:


# Criteria 3
# Min Free energy of pgRNA > -5
list3 = []
for val in list1:
    crRNA_combined = val
    list3.append(crRNA_combined)


# In[65]:


# Converting the list3 to a text file
with open('crRNA_combined.txt', 'w') as f:
    for item in list3:
        f.write("%s\n" % item)


# In[66]:


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


# In[67]:


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


# In[68]:


### Function to calculate the Secondary Structure Using RNAFold Software
def cal_Structure(s):
    s1 = s.split(" ")
    Structure = s1[0]
    return Structure


# In[69]:


command1 = "RNAfold < out_crRNA_combined.fsa >> energyoutput.txt"
subprocess.call(command1, shell=True) 


# In[70]:


command2 = "del *_ss.ps"
subprocess.call(command2, shell=True) 


# In[71]:


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


# In[72]:


energy_list = []
Structure_list = []
for val in energy_count:
    energy = cal_Min_Eng(val)
    energy_list.append(energy)
    Structure = cal_Structure(val)
    Structure_list.append(Structure)


# In[73]:


df_crRNA1["Free Energy"] = energy_list
df_crRNA1["Secondary Structure"] = Structure_list


# In[74]:


df_crRNA2 = df_crRNA1.loc[(df_crRNA1['Free Energy']) > -5]
len( df_crRNA2)


# In[75]:


df_crRNA2 = df_crRNA2.reset_index(drop=True)
df_crRNA2


# In[76]:


df_crRNA2.to_csv (r'df_crRNA_withEnergy_filtered.csv', index = False, header=True)


# In[77]:


list_filtered = df_crRNA2['pgRNA'].values.tolist()
  
len(list_filtered)


# In[78]:


reverse_complement_rna_list = []
for val in list_filtered:
    # Taking the reverse complement of pg-DNA sequence
    s = revcom(val)
    s = s.replace('T', 'U')
    reverse_complement_rna_list.append(s)
    
    


# In[79]:


len(reverse_complement_rna_list)


# In[80]:


# Converting the reverse_complement_rna_list into a text file
with open('pgRNA_filtered.txt', 'w') as f:
    for item in reverse_complement_rna_list:
        f.write("%s\n" % item)


# In[81]:


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

# In[82]:


command = 'blastn -query out_pgRNA_filtered.fsa -task blastn-short -db C:/Users/18064/Documents/Blast/db/Human_NCBI_DNAdb -outfmt "7 qacc sacc qstart qend sstart send" -out pgRNA_NCBI_blast.out'
subprocess.call(command, shell=True)


# In[83]:


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

# Making a list of blast hits in human genome for the query sequences
new_hit_list =[]
for val in hits_list:
    x = val.split("\n")[0].split(" ")[1]
    x = int(x)
    new_hit_list.append(x)


# In[84]:


blast_hits = pd.DataFrame({'GuideNumber': new_query_list,'Hits': new_hit_list})
blast_hits.head()


# In[85]:


len(blast_hits)


# In[86]:


df_crRNA2["HumanBlastHits"] = new_hit_list


# In[87]:


df_crRNA2.to_csv (r'df_crRNA_withhits.csv', index = False, header=True)


# In[88]:


f=open("pgRNA_NCBI_blast.out", "r")
contents =f.readlines()
Blast_list = []
for i in range(len(contents)):
    target = contents[i].find("#")
    if (target == -1):
        Blast_list.append(i)


# In[89]:


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


# In[90]:


# This is for Cas9
range_list = []
seq_4_range_list = []
l = len(sstart_list)
for i in range(l):
    if (sstart_list[i] > send_list[i]): # When hit top strand
        #print('strand = minus')

        range_blast = accession_list[i] + '  '+ str(int(sstart_list[i])+int(qstart_list[i])-1-23) + '-' + str(int(sstart_list[i])-1 + int(qend_list[i])) + '  '+ 'minus'
    else: # When hit Bottom strand
        range_blast = accession_list[i] + '  '+ str(int(sstart_list[i])+ 1 - int(qstart_list[i])) + '-' + str(int(sstart_list[i])-int(qstart_list[i])+1+23) + '  '+ 'plus'
    range_list.append(range_blast)
    seq_4_range_list.append(seq_list[i])


# In[91]:


range_list


# In[92]:


len(range_list)


# In[93]:


#seq_4_range_list


# In[94]:


# Converting the list to a text file
with open('range_NCBIBlasthits.txt', 'w') as f:
    for item in range_list:
        f.write("%s\n" % item)


# In[95]:


command1 = 'blastdbcmd -db C:/Users/18064/Documents/Blast/db/Human_NCBI_DNAdb -entry_batch range_NCBIBlasthits.txt -outfmt "%s" -out NCBI_Human_DNA_Range.txt'
subprocess.call(command1, shell=True)


# In[96]:


## Reading the NCBI_Human_RNA_Range.txt file
f_rangeDNA=open("NCBI_Human_DNA_Range.txt", "r")
contents_rangeDNA =f_rangeDNA.readlines()
Human_rangeDNA_list = []
for val in contents_rangeDNA:
    val1 = (val.split('\n')[0])[0:23]
    Human_rangeDNA_list.append(val1)   


# In[97]:


len(Human_rangeDNA_list)


# In[98]:


len(Human_rangeDNA_list[1])# Checking if there are 23 nucleotide in human DNA range


# In[99]:


seq_blast_cfd = []
GuideScore_blast = []
for val in seq_4_range_list:
    i = (int(val) - 1)
    #print(i)
    seq_blast = revcom(df_crRNA2.iloc[i].loc['pgRNA'])
    guide_score = df_crRNA2.iloc[i].loc['GuideScore']
    seq_blast_cfd.append(seq_blast)
    GuideScore_blast.append(guide_score)


# In[100]:


Blast_cfd_df = pd.DataFrame({'GuideNumber': seq_4_range_list,'SequenceComplementGuide': seq_blast_cfd,'HumanDNA_range': Human_rangeDNA_list,'GuideActivity': GuideScore_blast})


# In[101]:


Blast_cfd_df.head()
    


# #### Modified Function for CFD Mismatch score in Human Range, Now we have to consider the PAM effect

# In[102]:


# a is the Spacer gDNA/RNA and b is the pgDNA
def Cal_Mismatch_Score_blastcmd(a,b):
    wt = a[:-3] 
    off = b
    sg = off[:-3]
    pam = off[-2:] 
    score = calc_cfd(wt,sg,pam)
    return (score)


# In[103]:


## Calculating the cfd score for the Human transcriptome region corresponding to blast hits
x = len(Human_rangeDNA_list)
cfd_humanBlast_list = []
Activity_humanBlast_list = []
for i in range(x):
    #print(i)
    try:
        #cfd_humanBlast = Cal_Mismatch_Score(Blast_cfd_df["SequenceComplementGuide"][i], Blast_cfd_df["HumanDNA_range"][i])
        cfd_humanBlast = Cal_Mismatch_Score_blastcmd(Blast_cfd_df["HumanDNA_range"][i], Blast_cfd_df["SequenceComplementGuide"][i])
        Activity_humansite = cfd_humanBlast*Blast_cfd_df["GuideActivity"][i]
    except IndexError:
        cfd_humanBlast = 'null'
        Activity_humansite = "null"
    cfd_humanBlast_list.append(cfd_humanBlast)
    Activity_humanBlast_list.append(Activity_humansite)


# In[104]:


Cal_Mismatch_Score_blastcmd('ACACTTCCACCACTGTACTTTCT','CCACTTCCACCACTGTACTG')


# In[105]:


Blast_cfd_df["CFD_Score"] = cfd_humanBlast_list
Blast_cfd_df["ActivityScore"] = Activity_humanBlast_list


# In[106]:


Blast_cfd_df.head()


# In[107]:


Blast_cfd_df["Guide"] = Blast_cfd_df["SequenceComplementGuide"].apply(revcom)


# In[108]:


Blast_cfd_df.head()


# In[109]:


Blast_cfd_df.to_csv (r'BlastHumanRnaHits_Activity.csv', index = False, header=True)


# In[ ]:




