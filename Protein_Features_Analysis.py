#!/usr/bin/env python
# coding: utf-8

# In[12]:


from Bio import SeqIO
from Bio import SwissProt
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# In[1]:


#from all wild aa seq, filter all wild human aa seq
hum_aa_rec = []
i = 1

for seq_record in SeqIO.parse("uniprot_sprot.fasta", "fasta"):
    if '_HUMAN' in seq_record.id:
        record = SeqRecord(seq_record.seq, seq_record.id)
        hum_aa_rec.append(record)
        i += 1
        print(i)
        
#File with all wild human aa seq:
SeqIO.write(hum_aa_rec, "human_aa.fasta", "fasta")


# In[18]:


#Convert humsavar into a pd dataframe
import pandas as pd
from pandas import DataFrame
import re

data = pd.read_csv('filtered_humsavar.csv') #This is a csv file of all the entries in the humsavar database

#removing the 'p.' in AAC column
aac = data['AAC']
aac_f = []
for item in aac:
    aac_f.append(item[2:])
    
#change LB/B and LP/P to 0(benign) and 1(pathogenic)
cate = []
for row in data.index:
    if data.loc[row,'CATE'] == 'LB/B':
        cate.append(0)
    if data.loc[row,'CATE'] == 'LP/P':
        cate.append(1)
print(cate[:100])
    
#the dataframe for useful columns in humsavar:  
dict = {'gene':data['GENE'],'spac':data['SPAC'],'aac':aac_f,'cate':cate}
df = pd.DataFrame(dict)
df.to_csv('unshuffled_humsavar.csv')

df.head()


# In[22]:


#DO NOT RUN THIS STEP REPEATEDLY
#randomly shuffle the dataframe
df_shu = df.sample(frac=1)
df_shu.to_csv('shuffled_humsavar.csv')

df_shu.head()


# In[50]:


#filter out the spac of the first 20,000 benign and pathogenic mutations
i = 1
j = 1
wild_b_spac = [] #spac for first 20,000 benign mutations
wild_p_spac = [] #spac for first 20,000 pathogenic mutations
aac_b = [] #aa change for first 20,000 benign mutations
aac_p = [] #aa change for first 20,000 pathogenic mutations

df_shu = pd.read_csv('shuffled_humsavar.csv')

for row in df_shu.index:
    if df_shu.loc[row,'cate'] == 0 and i <= 20000:
        wild_b_spac.append(df_shu.loc[row,'spac'])
        aac_b.append(df_shu.loc[row,'aac'])
        i += 1
    if df_shu.loc[row,'cate'] == 1 and j <= 20000:    
        wild_p_spac.append(df_shu.loc[row,'spac'])
        aac_p.append(df_shu.loc[row,'aac'])
        j += 1

print(wild_b_spac[:5])
print(wild_p_spac[:5])

print(aac_b[:5])
print(aac_p[:5])

print(len(wild_b_spac)==len(aac_b)==20000)
print(len(wild_p_spac)==len(aac_p)==20000)


# In[98]:


#redo the step above with seqs less than 1022 residues KEY ASSUMPTION!!
rec_b_1022_f5 = []
rec_b_1022_l5 = []
rec_p_1022_f5 = []
rec_p_1022_l5 = []

#account for the indexes of which spac are used, to later select the corresponding aac
wb_spac_i = []#有一个重复，在[4998]和[4999]是一样的，因为应该是5000<i,j<=10000，但是不会影响最后的结果（一万条里面有一条重复），所以没改
wp_spac_i = []

i = 1
j = 1
for seq_record in SeqIO.parse("human_aa.fasta", "fasta"):
    full_id = seq_record.id.split('|') #full_id[1] corresponds to the spac format
    for index, item in enumerate(wild_b_spac):
        if full_id[1] == item and 10 < len(seq_record) < 1022 and i < 5000 and 'U' not in seq_record:
            rec_b_1022_f5.append(SeqRecord(seq_record.seq, seq_record.id))
            wb_spac_i.append(index)
            i += 1
        if full_id[1] == item and 10 < len(seq_record) < 1022 and 5000 <= i < 10000 and 'U' not in seq_record:
            rec_b_1022_l5.append(SeqRecord(seq_record.seq, seq_record.id))
            wb_spac_i.append(index)
            i += 1
                
    for index, item in enumerate(wild_p_spac):
        if full_id[1] == item and 10 < len(seq_record) < 1022 and j < 5000 and 'U' not in seq_record:
            rec_p_1022_f5.append(SeqRecord(seq_record.seq, seq_record.id))
            wp_spac_i.append(index)
            j += 1
        if full_id[1] == item and 10 < len(seq_record) < 1022 and 5000 <= j < 10000 and 'U' not in seq_record:
            rec_p_1022_l5.append(SeqRecord(seq_record.seq, seq_record.id))
            wp_spac_i.append(index)
            j += 1

print(len(wb_spac_i)==len(wp_spac_i)==9999)            
            
#writing these sequence in fasta file to input in NetSurfP-3.0
SeqIO.write(rec_b_1022_f5, "wild_aa_b_les1022res_f5.fasta", "fasta") #first 5000 wild benign amino acid sequences
SeqIO.write(rec_b_1022_l5, "wild_aa_b_les1022res_l5.fasta", "fasta") #last 5000 wild benign amino acid sequences
SeqIO.write(rec_p_1022_f5, "wild_aa_p_les1022res_f5.fasta", "fasta") #first 5000 wild pathogenic amino acid sequences
SeqIO.write(rec_p_1022_l5, "wild_aa_p_les1022res_l5.fasta", "fasta") #last 5000 wild pathogenic amino acid sequences


# In[97]:


#location of mutation in the benign and pathogenic sequences:
aac_b_used = [aac_b[x] for x in wb_spac_i]
aac_p_used = [aac_p[x] for x in wp_spac_i]

aac_int_b = [] #the int in each str in the aac column of humsavar
aac_int_p = []

for item in aac_b_used:
    aac_int_b.append([int(s) for s in re.findall(r'-?\d+\.?\d*', item)][0])

for item in aac_p_used:
    aac_int_p.append([int(s) for s in re.findall(r'-?\d+\.?\d*', item)][0])

print(aac_int_b[:5]) 
print(aac_int_p[:5])

#Finding the corresponding mutant aa seqs
rec_bm_1022_f5 = []
rec_bm_1022_l5 = []
rec_pm_1022_f5 = []
rec_pm_1022_l5 = []

i = 0
for seq_record in SeqIO.parse('wild_aa_b_les1022res_f5.fasta','fasta'):
    wil_seq = seq_record.seq 
    index = aac_int_b[i] - 1 #in python numbering
    mut_seq = wil_seq[:index] + aac_b_used[i][-3] + wil_seq[index+1:]
    rec_bm_1022_f5.append(SeqRecord(Seq(mut_seq),seq_record.id))
    i += 1

for seq_record in SeqIO.parse('wild_aa_b_les1022res_l5.fasta','fasta'):
    wil_seq = seq_record.seq 
    index = aac_int_b[i] - 1
    mut_seq = wil_seq[:index] + aac_b_used[i][-3] + wil_seq[index+1:]
    rec_bm_1022_l5.append(SeqRecord(Seq(mut_seq),seq_record.id))
    i += 1

j = 0
for seq_record in SeqIO.parse('wild_aa_p_les1022res_f5.fasta','fasta'):
    wil_seq = seq_record.seq 
    index = aac_int_p[j] - 1
    mut_seq = wil_seq[:index] + aac_p_used[j][-3] + wil_seq[index+1:]
    rec_pm_1022_f5.append(SeqRecord(Seq(mut_seq),seq_record.id))
    j += 1

for seq_record in SeqIO.parse('wild_aa_p_les1022res_l5.fasta','fasta'):
    wil_seq = seq_record.seq 
    index = aac_int_p[j] - 1
    mut_seq = wil_seq[:index] + aac_p_used[j][-3] + wil_seq[index+1:]
    rec_pm_1022_l5.append(SeqRecord(Seq(mut_seq),seq_record.id))
    j += 1
    
#writing these sequence in fasta file to input in NetSurfP-3.0
SeqIO.write(rec_bm_1022_f5, "mut_aa_b_les1022res_f5.fasta", "fasta") #first 5000 wild benign amino acid sequences
SeqIO.write(rec_bm_1022_l5, "mut_aa_b_les1022res_l5.fasta", "fasta") #last 5000 wild benign amino acid sequences
SeqIO.write(rec_pm_1022_f5, "mut_aa_p_les1022res_f5.fasta", "fasta") #first 5000 wild pathogenic amino acid sequences
SeqIO.write(rec_pm_1022_l5, "mut_aa_p_les1022res_l5.fasta", "fasta") #last 5000 wild pathogenic amino acid sequences


# In[101]:


#check the number of overlaps in the wild fasta files
overlap = []
for item in wild_b_spac:
    if item in wild_p_spac:
        overlap.append(item)
        
print('Out of 40,000 wilds sequences, ' + str(len(overlap)*2) + ' of them overlapped, which is ' + 
      str(len(overlap)*2/40000*100) + '% of the sequences.')
#This overlap will only affect the wild (benign and pathongeic) sequences, which will be taken care by considering delta_rsa

#now, I can feed the sequences into NetsurfP-3.0.


# In[43]:


#By now, I have obtained the .csv files from NetsurfP-3.0
#Since aa seq only accounts for first digit, there must be cases where the mutation is unapparent.

import pandas as pd
from pandas import DataFrame

#creating a pd dataframe for the benign mutations
wb1 = pd.read_csv('wb_f5.csv')
wb2 = pd.read_csv('wb_l5.csv')

wb_t = wb1.append(wb2)

ID = wb_t['id']
SEQ_w = wb_t[' seq']
N = wb_t[' n']
RSA_w = wb_t[' rsa']

mb1 = pd.read_csv('mb_f5.csv')
mb2 = pd.read_csv('mb_l5.csv')

mb_t = mb1.append(mb2)

SEQ_m = mb_t[' seq']
RSA_m = mb_t[' rsa']
Delta_RSA = RSA_w - RSA_m

dict = {'ID':ID, 'N':N, 'J':N, 'seq_w':SEQ_w, 'seq_m':SEQ_m, 'rsa_w':RSA_w, 'rsa_m':RSA_m, 'delta_rsa':Delta_RSA}
df = pd.DataFrame(dict)

df.to_csv('benign_without_j.csv')
print(len(wb1)+len(wb2)==len(wb_t))

print(df.head())


# In[44]:


#removing the repeated mutations -- room for error (e.g. removed the correct one and kept the wrong one)
df = pd.read_csv('benign_without_j.csv')

t2 = [] #index of the first residue in every sequence
for row in df.index:
    if df.loc[row,'N'] == 1:
        t2.append(row)

#identifying all the repeated sequences at the location of mutation 
i = 1
k = 0
t_2 = 0
rep_rows = []
for row in df.index:
    if i <= t2[t_2 + 1] and df.loc[row, 'seq_m'] == df.loc[row, 'seq_w']:
        k += 1
        if k == t2[t_2 + 1] - t2[t_2]:
            rep_rows.append(row)
    if i == t2[t_2 + 1] and t_2 < len(t2) - 2:
        t_2 += 1
        k = 0
    i += 1
    
print('the number of repeated rows is ' + str(len(rep_rows)))
print(rep_rows[:5])


# In[45]:


#finding the first and last residues of these repeated sequences
ALL1 = [] #last residue
ALL2 = [] #first residue
for item in rep_rows:
    run = False
    for i in range(len(t2)):
        if t2[i] >= item and not run:
            ALL1.append(t2[i])
            ALL2.append(t2[i-1])
            run = True

#removing all these residues            
i = 0
for item in ALL1:
    DROP = list(range(ALL2[i], ALL1[i]))
    df.drop(DROP,axis = 0,inplace = True)
    i += 1
    
df2 = df.reset_index()
print('deleted all repeated sequences')

#adding j values
t1 = [] #a list of accumulated rows where J=0
t2 = [] #a list of accumulated rows where n=1
for row in df2.index:
    w_seq = df2.loc[row, 'seq_w']   
    m_seq = df2.loc[row, 'seq_m']
    nn = df2.loc[row, 'N']
    if w_seq != m_seq:
        df2.loc[row,'J'] = 0
        t1.append(row)
    if w_seq == m_seq and nn != 1:
        df2.loc[row,'J'] = 'NA'
    if nn == 1 and row != 0:
        df2.loc[row,'J'] = 'MAX'
        t2.append(row)

t_1 = 0
t_2 = 0
i = 0
J4 = [] #column of j values
for row in df2.index:
    J4.append( t1[t_1] - i )
    if i == t2[t_2] - 1 and t_2 < len(t2) - 1:
        t_1 = t_1 + 1
        t_2 = t_2 + 1
    i = i + 1         

df2['J'] = J4
df2.to_csv('benign_filtered_data.csv')
print('csv file created')


# In[49]:


#repeat for pathogenic sequences
import pandas as pd
from pandas import DataFrame

#creating a pd dataframe for the benign mutations
wp1 = pd.read_csv('wp_f5.csv')
wp2 = pd.read_csv('wp_l5.csv')

wp_t = wp1.append(wp2)

ID = wp_t['id']
SEQ_w = wp_t[' seq']
N = wp_t[' n']
RSA_w = wp_t[' rsa']

mp1 = pd.read_csv('mp_f5.csv')
mp2 = pd.read_csv('mp_l5.csv')

mp_t = mp1.append(mp2)

SEQ_m = mp_t[' seq']
RSA_m = mp_t[' rsa']
Delta_RSA = RSA_w - RSA_m

dict = {'ID':ID, 'N':N, 'J':N, 'seq_w':SEQ_w, 'seq_m':SEQ_m, 'rsa_w':RSA_w, 'rsa_m':RSA_m, 'delta_rsa':Delta_RSA}
df = pd.DataFrame(dict)

print(len(wb1)+len(wb2)==len(wb_t))
df.to_csv('pathogenic_without_j.csv')
df.head()


# In[50]:


df = pd.read_csv('pathogenic_without_j.csv')

t2 = [] #index of the first residue in every sequence
for row in df.index:
    if df.loc[row,'N'] == 1:
        t2.append(row)

#identifying all the repeated sequences at the location of mutation 
i = 1
k = 0
t_2 = 0
rep_rows = []
for row in df.index:
    if i <= t2[t_2 + 1] and df.loc[row, 'seq_m'] == df.loc[row, 'seq_w']:
        k += 1
        if k == t2[t_2 + 1] - t2[t_2]:
            rep_rows.append(row)
    if i == t2[t_2 + 1] and t_2 < len(t2) - 2:
        t_2 += 1
        k = 0
    i += 1
    
print('the number of repeated rows is ' + str(len(rep_rows)))
print(rep_rows[:5])

#finding the first and last residues of these repeated sequences
ALL1 = [] #last residue
ALL2 = [] #first residue
for item in rep_rows:
    run = False
    for i in range(len(t2)):
        if t2[i] >= item and not run:
            ALL1.append(t2[i])
            ALL2.append(t2[i-1])
            run = True

#removing all these residues            
i = 0
for item in ALL1:
    DROP = list(range(ALL2[i], ALL1[i]))
    df.drop(DROP,axis = 0,inplace = True)
    i += 1
    
df2 = df.reset_index()
print('deleted all repeated sequences')

#adding j values
t1 = [] #a list of accumulated rows where J=0
t2 = [] #a list of accumulated rows where n=1
for row in df2.index:
    w_seq = df2.loc[row, 'seq_w']   
    m_seq = df2.loc[row, 'seq_m']
    nn = df2.loc[row, 'N']
    if w_seq != m_seq:
        df2.loc[row,'J'] = 0
        t1.append(row)
    if w_seq == m_seq and nn != 1:
        df2.loc[row,'J'] = 'NA'
    if nn == 1 and row != 0:
        df2.loc[row,'J'] = 'MAX'
        t2.append(row)

t_1 = 0
t_2 = 0
i = 0
J4 = [] #column of j values
for row in df2.index:
    J4.append( t1[t_1] - i )
    if i == t2[t_2] - 1 and t_2 < len(t2) - 1:
        t_1 = t_1 + 1
        t_2 = t_2 + 1
    i = i + 1         

df2['J'] = J4
df2.to_csv('pathogenic_filtered_data.csv')
print('csv file created')


# In[6]:


#Computing the initial rsa data at variant location
import pandas as pd
from pandas import DataFrame

bb = pd.read_csv('benign_filtered_data.csv')
pp = pd.read_csv('pathogenic_filtered_data.csv')

rsa_i_b = []
for row in bb.index:
    if bb.loc[row,'seq_w'] != bb.loc[row,'seq_m']:
        rsa_i_b.append(bb.loc[row, 'rsa_w'])


rsa_i_p = []
for row in pp.index:
    if pp.loc[row,'seq_w'] != pp.loc[row,'seq_m']:
        rsa_i_p.append(pp.loc[row, 'rsa_w'])


# In[7]:


#Graph 1: distribution of initial rsa at the variant location
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.distplot(rsa_i_b, hist = True, kde = False, bins = 40, 
                 kde_kws = {'linewidth': 2}, color = 'b')
sns.distplot(rsa_i_p, hist = True, kde = False, bins = 40, 
                 kde_kws = {'linewidth': 2}, color = 'r')

plt.title("Distribution of initial rsa values at variant site")
plt.xlabel("Initial rsa values")
plt.ylabel("Frequency")
plt.legend(["benign", "pathogenic"], loc=0)
plt.xlim(0,1)


# In[8]:


#Computing the final rsa data at variant location
import pandas as pd
from pandas import DataFrame

bb = pd.read_csv('benign_filtered_data.csv')
pp = pd.read_csv('pathogenic_filtered_data.csv')

rsa_f_b = []
for row in bb.index:
    if bb.loc[row,'seq_w'] != bb.loc[row,'seq_m']:
        rsa_f_b.append(bb.loc[row, 'rsa_m'])


rsa_f_p = []
for row in pp.index:
    if pp.loc[row,'seq_w'] != pp.loc[row,'seq_m']:
        rsa_f_p.append(pp.loc[row, 'rsa_m'])


# In[9]:


#Graph 2: distribution of final rsa at the variant location
sns.distplot(rsa_f_b, hist = True, kde = False, bins = 40, 
                 kde_kws = {'linewidth': 2}, color = 'b')
sns.distplot(rsa_f_p, hist = True, kde = False, bins = 40, 
                 kde_kws = {'linewidth': 2}, color = 'r')

plt.title("Distribution of final rsa values at variant site")
plt.xlabel("Final rsa values")
plt.ylabel("Frequency")
plt.legend(["benign", "pathogenic"], loc=0)
plt.xlim(0,1)


# In[20]:


#Graph 3: J distribution graph in both directions
import pandas as pd
import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt

x_range = range(min(bb['J']), max(bb['J']+1))
x = np.array(list(x_range)) # a list from -991 to 1014

print(max(bb['J']))
print(min(bb['J']))
print(max(pp['J']))
print(min(pp['J']))

J_P = [] #a list of the missing, positive J values adding to the pathogenic data because it has shorter residues
for i in range(max(pp['J']),max(bb['J'])+1):
    J_P.append(i)
    
J_N = [] #the corresponding list of the missing, negative J values
for i in range(abs(min(pp['J'])),abs(min(bb['J']))+1):
    J_N.append(-i)
J_add = J_P + J_N    
    
NA_rsa = []
for i in range(len(J_add)):
    NA_rsa.append(None) 

pp2 = pd.DataFrame({'J':pp['J'],'delta_rsa':pp['delta_rsa']})

dict1 = {'J':J_add, 'delta_rsa':NA_rsa}
df1 = pd.DataFrame(dict1)
                
pp3 = pp2.append(df1) #creating a supplemented pathogenic dataset
f_p = pp3.groupby('J').mean()

bb2 = pd.DataFrame({'J':bb['J'],'delta_rsa':bb['delta_rsa']})
f_b = bb2.groupby('J').mean()

y = np.array(f_b['delta_rsa'])
z = np.array(f_p['delta_rsa'])

#the asbolute delta rsa values:
y_abs = np.array(abs(f_b['delta_rsa'])) 
z_abs = np.array(abs(f_p['delta_rsa']))

plt.plot(x,y_abs,'r')
plt.plot(x,z_abs,'b')
plt.title('Distribution of absolute delta_rsa values at every value of J')
plt.xlabel('J values')
plt.ylabel('Absolute delta_rsa values')
plt.xlim(-1100,1100)
plt.legend(["benign", "pathogenic"], loc=0)
plt.show()

plt.plot(x,y_abs,'r')
plt.plot(x,z_abs,'b')
plt.title('Distribution of absolute delta_rsa values at every value of J')
plt.xlabel('J values')
plt.ylabel('Absolute delta_rsa values')
plt.xlim(-75,75)
plt.legend(["benign", "pathogenic"], loc=0)
plt.show()


# In[10]:


#Computing delta_rsa data
import pandas as pd
from pandas import DataFrame

bb = pd.read_csv('benign_filtered_data.csv')
pp = pd.read_csv('pathogenic_filtered_data.csv')

d_rsa_b = []
for row in bb.index:
    if bb.loc[row,'seq_w'] != bb.loc[row,'seq_m']:
        d_rsa_b.append(bb.loc[row, 'delta_rsa'])


d_rsa_p = []
for row in pp.index:
    if pp.loc[row,'seq_w'] != pp.loc[row,'seq_m']:
        d_rsa_p.append(pp.loc[row, 'delta_rsa'])


# In[14]:


#Graph 4: Distribution of the average delta_rsa value at 20 residues near by the variant site
bb = pd.read_csv('benign_filtered_data.csv')
pp = pd.read_csv('pathogenic_filtered_data.csv')

x_1 = '(-10^-1,-10^-2]'
x_2 = '(-10^-2,-10^-3]'
x_3 = '(-10^-3,-10^-4]'
x_4 = '(-10^-4,-10^-5]'
x_5 = '(-10^-5,-10^-6]'
x_6 = '(-10^-6,0]'
x_7 = '(0,10^-6]'
x_8 = '(10^-6,10^-5]'
x_9 = '(10^-5,10^-4]'
x_10 = '(10^-4,10^-3]'
x_11 = '(10^-3,10^-2]'
x_12 = '(10^-2,10^-1]'

x = np.array([x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9,x_10,x_11,x_12])

J_AVE_b = []
for row in bb.index:
    if bb.loc[row,'seq_w'] != bb.loc[row,'seq_m']:
        J_AVE_b.append((bb.loc[row, 'delta_rsa']+bb.loc[row+1, 'delta_rsa']+bb.loc[row+2, 'delta_rsa']
                     +bb.loc[row+3, 'delta_rsa']+bb.loc[row+4, 'delta_rsa']+bb.loc[row+5, 'delta_rsa']
                     +bb.loc[row+6, 'delta_rsa']+bb.loc[row+7, 'delta_rsa']+bb.loc[row+8, 'delta_rsa']
                     +bb.loc[row+9, 'delta_rsa']+bb.loc[row-1, 'delta_rsa']+bb.loc[row-2, 'delta_rsa']
                     +bb.loc[row-3, 'delta_rsa']+bb.loc[row-4, 'delta_rsa']+bb.loc[row-5, 'delta_rsa']
                     +bb.loc[row-6, 'delta_rsa']+bb.loc[row-7, 'delta_rsa']+bb.loc[row-8, 'delta_rsa']
                     +bb.loc[row-9, 'delta_rsa'])/19)

bin1 = []
bin2 = []
bin3 = []
bin4 = []
bin5 = []
bin6 = []
bin7 = []
bin8 = []
bin9 = []
bin10 = []
bin11 = []
bin12 = []
for i in J_AVE_b:
    if -10**-1 < i <= -10**-2:
        bin1.append(i)
    if -10**-2 < i <= -10**-3:
        bin2.append(i)
    if -10**-3 < i <= -10**-4:
        bin3.append(i)
    if -10**-4 < i <= -10**-5:
        bin4.append(i)
    if -10**-5 < i <= -10**-6:
        bin5.append(i)
    if -10**-6 < i <= 0:
        bin6.append(i)
    if 0 < i <= 10**-6:
        bin7.append(i)
    if 10**-6 < i <= 10**-5:
        bin8.append(i)
    if 10**-5 < i <= 10**-4:
        bin9.append(i)
    if 10**-4 < i <= 10**-3:
        bin10.append(i)
    if 10**-3 < i <= 10**-2:
        bin11.append(i)
    if 10**-2 < i <= 10**-1:
        bin12.append(i)
        
# print(len(bin1)+len(bin2)+len(bin3)+len(bin4)+len(bin5)+len(bin6)+len(bin7)+len(bin8)+len(bin9)+len(bin10)
#       +len(bin11)+len(bin12)==len(J_AVE_b))

y_1 = len(bin1)/len(J_AVE_b)
y_2 = len(bin2)/len(J_AVE_b)
y_3 = len(bin3)/len(J_AVE_b)
y_4 = len(bin4)/len(J_AVE_b)
y_5 = len(bin5)/len(J_AVE_b)
y_6 = len(bin6)/len(J_AVE_b)
y_7 = len(bin7)/len(J_AVE_b)
y_8 = len(bin8)/len(J_AVE_b)
y_9 = len(bin9)/len(J_AVE_b)
y_10 = len(bin10)/len(J_AVE_b)
y_11 = len(bin11)/len(J_AVE_b)
y_12 = len(bin12)/len(J_AVE_b)

y = np.array([y_1,y_2,y_3,y_4,y_5,y_6,y_7,y_8,y_9,y_10,y_11,y_12])


J_AVE_p = []
for row in pp.index:
    if pp.loc[row,'seq_w'] != pp.loc[row,'seq_m']:
        J_AVE_p.append((pp.loc[row, 'delta_rsa']+pp.loc[row+1, 'delta_rsa']+pp.loc[row+2, 'delta_rsa']
                     +pp.loc[row+3, 'delta_rsa']+pp.loc[row+4, 'delta_rsa']+pp.loc[row+5, 'delta_rsa']
                     +pp.loc[row+6, 'delta_rsa']+pp.loc[row+7, 'delta_rsa']+pp.loc[row+8, 'delta_rsa']
                     +pp.loc[row+9, 'delta_rsa']+pp.loc[row-1, 'delta_rsa']+pp.loc[row-2, 'delta_rsa']
                     +pp.loc[row-3, 'delta_rsa']+pp.loc[row-4, 'delta_rsa']+pp.loc[row-5, 'delta_rsa']
                     +pp.loc[row-6, 'delta_rsa']+pp.loc[row-7, 'delta_rsa']+pp.loc[row-8, 'delta_rsa']
                     +pp.loc[row-9, 'delta_rsa'])/19)

bin1 = []
bin2 = []
bin3 = []
bin4 = []
bin5 = []
bin6 = []
bin7 = []
bin8 = []
bin9 = []
bin10 = []
bin11 = []
bin12 = []
for i in J_AVE_p:
    if -10**-1 < i <= -10**-2:
        bin1.append(i)
    if -10**-2 < i <= -10**-3:
        bin2.append(i)
    if -10**-3 < i <= -10**-4:
        bin3.append(i)
    if -10**-4 < i <= -10**-5:
        bin4.append(i)
    if -10**-5 < i <= -10**-6:
        bin5.append(i)
    if -10**-6 < i <= 0:
        bin6.append(i)
    if 0 < i <= 10**-6:
        bin7.append(i)
    if 10**-6 < i <= 10**-5:
        bin8.append(i)
    if 10**-5 < i <= 10**-4:
        bin9.append(i)
    if 10**-4 < i <= 10**-3:
        bin10.append(i)
    if 10**-3 < i <= 10**-2:
        bin11.append(i)
    if 10**-2 < i <= 10**-1:
        bin12.append(i)
        
# print(len(bin1)+len(bin2)+len(bin3)+len(bin4)+len(bin5)+len(bin6)+len(bin7)+len(bin8)+len(bin9)+len(bin10)
#       +len(bin11)+len(bin12)==len(J_AVE_p))

z_1 = len(bin1)/len(J_AVE_p)
z_2 = len(bin2)/len(J_AVE_p)
z_3 = len(bin3)/len(J_AVE_p)
z_4 = len(bin4)/len(J_AVE_p)
z_5 = len(bin5)/len(J_AVE_p)
z_6 = len(bin6)/len(J_AVE_p)
z_7 = len(bin7)/len(J_AVE_p)
z_8 = len(bin8)/len(J_AVE_p)
z_9 = len(bin9)/len(J_AVE_p)
z_10 = len(bin10)/len(J_AVE_p)
z_11 = len(bin11)/len(J_AVE_p)
z_12 = len(bin12)/len(J_AVE_p)

z = np.array([z_1,z_2,z_3,z_4,z_5,z_6,z_7,z_8,z_9,z_10,z_11,z_12])

fig, ax = plt.subplots()
ax.plot(x,y,'r') #benign
ax.plot(x,z,'b') #pathogenic
plt.title("probability distribution of average delta_rsa values")
plt.xlabel("X = average delta_rsa")
plt.ylabel("P(X)")
plt.legend(["benign", "pathogenic"], loc=0)

fig.autofmt_xdate()
plt.show()


# In[4]:


#Graph 5: Distribution of the relative average delta_rsa value at 20 residues near by the variant site
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from statistics import mean


bb = pd.read_csv('benign_filtered_data.csv')
pp = pd.read_csv('pathogenic_filtered_data.csv')

x_1 = '(-10^-1,-10^-2]'
x_2 = '(-10^-2,-10^-3]'
x_3 = '(-10^-3,-10^-4]'
x_4 = '(-10^-4,-10^-5]'
x_5 = '(-10^-5,-10^-6]'
x_6 = '(-10^-6,0]'
x_7 = '(0,10^-6]'
x_8 = '(10^-6,10^-5]'
x_9 = '(10^-5,10^-4]'
x_10 = '(10^-4,10^-3]'
x_11 = '(10^-3,10^-2]'
x_12 = '(10^-2,10^-1]'

x = np.array([x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9,x_10,x_11,x_12])

J_AVE_b = []
for row in bb.index:
    if bb.loc[row,'seq_w'] != bb.loc[row,'seq_m']:
        J_AVE_b.append(((bb.loc[row, 'delta_rsa']+bb.loc[row+1, 'delta_rsa']+bb.loc[row+2, 'delta_rsa']
                     +bb.loc[row+3, 'delta_rsa']+bb.loc[row+4, 'delta_rsa']+bb.loc[row+5, 'delta_rsa']
                     +bb.loc[row+6, 'delta_rsa']+bb.loc[row+7, 'delta_rsa']+bb.loc[row+8, 'delta_rsa']
                     +bb.loc[row+9, 'delta_rsa']+bb.loc[row-1, 'delta_rsa']+bb.loc[row-2, 'delta_rsa']
                     +bb.loc[row-3, 'delta_rsa']+bb.loc[row-4, 'delta_rsa']+bb.loc[row-5, 'delta_rsa']
                     +bb.loc[row-6, 'delta_rsa']+bb.loc[row-7, 'delta_rsa']+bb.loc[row-8, 'delta_rsa']
                     +bb.loc[row-9, 'delta_rsa'])/19)/((bb.loc[row, 'rsa_w']+bb.loc[row+1, 'rsa_w']+bb.loc[row+2, 'rsa_w']
                     +bb.loc[row+3, 'rsa_w']+bb.loc[row+4, 'rsa_w']+bb.loc[row+5, 'rsa_w']
                     +bb.loc[row+6, 'rsa_w']+bb.loc[row+7, 'rsa_w']+bb.loc[row+8, 'rsa_w']
                     +bb.loc[row+9, 'rsa_w']+bb.loc[row-1, 'rsa_w']+bb.loc[row-2, 'rsa_w']
                     +bb.loc[row-3, 'rsa_w']+bb.loc[row-4, 'rsa_w']+bb.loc[row-5, 'rsa_w']
                     +bb.loc[row-6, 'rsa_w']+bb.loc[row-7, 'rsa_w']+bb.loc[row-8, 'rsa_w']
                     +bb.loc[row-9, 'rsa_w'])))

bin1 = []
bin2 = []
bin3 = []
bin4 = []
bin5 = []
bin6 = []
bin7 = []
bin8 = []
bin9 = []
bin10 = []
bin11 = []
bin12 = []
for i in J_AVE_b:
    if -10**-1 < i <= -10**-2:
        bin1.append(i)
    if -10**-2 < i <= -10**-3:
        bin2.append(i)
    if -10**-3 < i <= -10**-4:
        bin3.append(i)
    if -10**-4 < i <= -10**-5:
        bin4.append(i)
    if -10**-5 < i <= -10**-6:
        bin5.append(i)
    if -10**-6 < i <= 0:
        bin6.append(i)
    if 0 < i <= 10**-6:
        bin7.append(i)
    if 10**-6 < i <= 10**-5:
        bin8.append(i)
    if 10**-5 < i <= 10**-4:
        bin9.append(i)
    if 10**-4 < i <= 10**-3:
        bin10.append(i)
    if 10**-3 < i <= 10**-2:
        bin11.append(i)
    if 10**-2 < i <= 10**-1:
        bin12.append(i)
        
# print(len(bin1)+len(bin2)+len(bin3)+len(bin4)+len(bin5)+len(bin6)+len(bin7)+len(bin8)+len(bin9)+len(bin10)
#       +len(bin11)+len(bin12)==len(J_AVE_b))

y_1 = len(bin1)/len(J_AVE_b)
y_2 = len(bin2)/len(J_AVE_b)
y_3 = len(bin3)/len(J_AVE_b)
y_4 = len(bin4)/len(J_AVE_b)
y_5 = len(bin5)/len(J_AVE_b)
y_6 = len(bin6)/len(J_AVE_b)
y_7 = len(bin7)/len(J_AVE_b)
y_8 = len(bin8)/len(J_AVE_b)
y_9 = len(bin9)/len(J_AVE_b)
y_10 = len(bin10)/len(J_AVE_b)
y_11 = len(bin11)/len(J_AVE_b)
y_12 = len(bin12)/len(J_AVE_b)

y = np.array([y_1,y_2,y_3,y_4,y_5,y_6,y_7,y_8,y_9,y_10,y_11,y_12])


J_AVE_p = []
for row in pp.index:
    if pp.loc[row,'seq_w'] != pp.loc[row,'seq_m']:
        J_AVE_p.append(((pp.loc[row, 'delta_rsa']+pp.loc[row+1, 'delta_rsa']+pp.loc[row+2, 'delta_rsa']
                     +pp.loc[row+3, 'delta_rsa']+pp.loc[row+4, 'delta_rsa']+pp.loc[row+5, 'delta_rsa']
                     +pp.loc[row+6, 'delta_rsa']+pp.loc[row+7, 'delta_rsa']+pp.loc[row+8, 'delta_rsa']
                     +pp.loc[row+9, 'delta_rsa']+pp.loc[row-1, 'delta_rsa']+pp.loc[row-2, 'delta_rsa']
                     +pp.loc[row-3, 'delta_rsa']+pp.loc[row-4, 'delta_rsa']+pp.loc[row-5, 'delta_rsa']
                     +pp.loc[row-6, 'delta_rsa']+pp.loc[row-7, 'delta_rsa']+pp.loc[row-8, 'delta_rsa']
                     +pp.loc[row-9, 'delta_rsa'])/19)/((pp.loc[row, 'rsa_w']+pp.loc[row+1, 'rsa_w']+pp.loc[row+2, 'rsa_w']
                     +pp.loc[row+3, 'rsa_w']+pp.loc[row+4, 'rsa_w']+pp.loc[row+5, 'rsa_w']
                     +pp.loc[row+6, 'rsa_w']+pp.loc[row+7, 'rsa_w']+pp.loc[row+8, 'rsa_w']
                     +pp.loc[row+9, 'rsa_w']+pp.loc[row-1, 'rsa_w']+pp.loc[row-2, 'rsa_w']
                     +pp.loc[row-3, 'rsa_w']+pp.loc[row-4, 'rsa_w']+pp.loc[row-5, 'rsa_w']
                     +pp.loc[row-6, 'rsa_w']+pp.loc[row-7, 'rsa_w']+pp.loc[row-8, 'rsa_w']
                     +pp.loc[row-9, 'rsa_w'])/19))


bin1 = []
bin2 = []
bin3 = []
bin4 = []
bin5 = []
bin6 = []
bin7 = []
bin8 = []
bin9 = []
bin10 = []
bin11 = []
bin12 = []
for i in J_AVE_p:
    if -10**-1 < i <= -10**-2:
        bin1.append(i)
    if -10**-2 < i <= -10**-3:
        bin2.append(i)
    if -10**-3 < i <= -10**-4:
        bin3.append(i)
    if -10**-4 < i <= -10**-5:
        bin4.append(i)
    if -10**-5 < i <= -10**-6:
        bin5.append(i)
    if -10**-6 < i <= 0:
        bin6.append(i)
    if 0 < i <= 10**-6:
        bin7.append(i)
    if 10**-6 < i <= 10**-5:
        bin8.append(i)
    if 10**-5 < i <= 10**-4:
        bin9.append(i)
    if 10**-4 < i <= 10**-3:
        bin10.append(i)
    if 10**-3 < i <= 10**-2:
        bin11.append(i)
    if 10**-2 < i <= 10**-1:
        bin12.append(i)
        
# print(len(bin1)+len(bin2)+len(bin3)+len(bin4)+len(bin5)+len(bin6)+len(bin7)+len(bin8)+len(bin9)+len(bin10)
#       +len(bin11)+len(bin12)==len(J_AVE_p))

z_1 = len(bin1)/len(J_AVE_p)
z_2 = len(bin2)/len(J_AVE_p)
z_3 = len(bin3)/len(J_AVE_p)
z_4 = len(bin4)/len(J_AVE_p)
z_5 = len(bin5)/len(J_AVE_p)
z_6 = len(bin6)/len(J_AVE_p)
z_7 = len(bin7)/len(J_AVE_p)
z_8 = len(bin8)/len(J_AVE_p)
z_9 = len(bin9)/len(J_AVE_p)
z_10 = len(bin10)/len(J_AVE_p)
z_11 = len(bin11)/len(J_AVE_p)
z_12 = len(bin12)/len(J_AVE_p)

z = np.array([z_1,z_2,z_3,z_4,z_5,z_6,z_7,z_8,z_9,z_10,z_11,z_12])

fig, ax = plt.subplots()
ax.plot(x,y,'r') #benign
ax.plot(x,z,'b') #pathogenic
plt.title("probability distribution of relative average delta_rsa values")
plt.xlabel("X = relative average delta_rsa")
plt.ylabel("P(X)")
plt.legend(["benign", "pathogenic"], loc=0)

fig.autofmt_xdate()
plt.show()


# In[25]:


#computing average delta rsa values for graph 6
J_AVE_b = []
for row in bb.index:
    if bb.loc[row,'seq_w'] != bb.loc[row,'seq_m']:
        J_AVE_b.append((bb.loc[row, 'delta_rsa']+bb.loc[row+1, 'delta_rsa']+bb.loc[row+2, 'delta_rsa']
                     +bb.loc[row+3, 'delta_rsa']+bb.loc[row+4, 'delta_rsa']+bb.loc[row+5, 'delta_rsa']
                     +bb.loc[row+6, 'delta_rsa']+bb.loc[row+7, 'delta_rsa']+bb.loc[row+8, 'delta_rsa']
                     +bb.loc[row+9, 'delta_rsa']+bb.loc[row-1, 'delta_rsa']+bb.loc[row-2, 'delta_rsa']
                     +bb.loc[row-3, 'delta_rsa']+bb.loc[row-4, 'delta_rsa']+bb.loc[row-5, 'delta_rsa']
                     +bb.loc[row-6, 'delta_rsa']+bb.loc[row-7, 'delta_rsa']+bb.loc[row-8, 'delta_rsa']
                     +bb.loc[row-9, 'delta_rsa'])/19)
        
J_AVE_p = []
for row in pp.index:
    if pp.loc[row,'seq_w'] != pp.loc[row,'seq_m']:
        J_AVE_p.append((pp.loc[row, 'delta_rsa']+pp.loc[row+1, 'delta_rsa']+pp.loc[row+2, 'delta_rsa']
                     +pp.loc[row+3, 'delta_rsa']+pp.loc[row+4, 'delta_rsa']+pp.loc[row+5, 'delta_rsa']
                     +pp.loc[row+6, 'delta_rsa']+pp.loc[row+7, 'delta_rsa']+pp.loc[row+8, 'delta_rsa']
                     +pp.loc[row+9, 'delta_rsa']+pp.loc[row-1, 'delta_rsa']+pp.loc[row-2, 'delta_rsa']
                     +pp.loc[row-3, 'delta_rsa']+pp.loc[row-4, 'delta_rsa']+pp.loc[row-5, 'delta_rsa']
                     +pp.loc[row-6, 'delta_rsa']+pp.loc[row-7, 'delta_rsa']+pp.loc[row-8, 'delta_rsa']
                     +pp.loc[row-9, 'delta_rsa'])/19)


# In[26]:


#Graph 6: Distribution of the average delta_rsa value at 20 residues near by the variant site - normal x bins
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.distplot(J_AVE_b, hist = True, kde = True, bins = 80, 
                 kde_kws = {'linewidth': 2}, color = 'b')
sns.distplot(J_AVE_p, hist = True, kde = True, bins = 80, 
                 kde_kws = {'linewidth': 2}, color = 'r')

plt.title("probability distribution of average delta_rsa values")
plt.xlabel("X = average delta_rsa")
plt.ylabel("P(x)")
plt.legend(["benign", "pathogenic"], loc=0)
plt.xlim(-0.05,0.05)


# In[21]:


#Graph 7: relative delta RSA for evenly distributed x bins
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from statistics import mean


bb = pd.read_csv('benign_filtered_data.csv')
pp = pd.read_csv('pathogenic_filtered_data.csv')

J_AVE_b = []
for row in bb.index:
    if bb.loc[row,'seq_w'] != bb.loc[row,'seq_m']:
        J_AVE_b.append(((bb.loc[row, 'delta_rsa']+bb.loc[row+1, 'delta_rsa']+bb.loc[row+2, 'delta_rsa']
                     +bb.loc[row+3, 'delta_rsa']+bb.loc[row+4, 'delta_rsa']+bb.loc[row+5, 'delta_rsa']
                     +bb.loc[row+6, 'delta_rsa']+bb.loc[row+7, 'delta_rsa']+bb.loc[row+8, 'delta_rsa']
                     +bb.loc[row+9, 'delta_rsa']+bb.loc[row-1, 'delta_rsa']+bb.loc[row-2, 'delta_rsa']
                     +bb.loc[row-3, 'delta_rsa']+bb.loc[row-4, 'delta_rsa']+bb.loc[row-5, 'delta_rsa']
                     +bb.loc[row-6, 'delta_rsa']+bb.loc[row-7, 'delta_rsa']+bb.loc[row-8, 'delta_rsa']
                     +bb.loc[row-9, 'delta_rsa'])/19)/((bb.loc[row, 'rsa_w']+bb.loc[row+1, 'rsa_w']+bb.loc[row+2, 'rsa_w']
                     +bb.loc[row+3, 'rsa_w']+bb.loc[row+4, 'rsa_w']+bb.loc[row+5, 'rsa_w']
                     +bb.loc[row+6, 'rsa_w']+bb.loc[row+7, 'rsa_w']+bb.loc[row+8, 'rsa_w']
                     +bb.loc[row+9, 'rsa_w']+bb.loc[row-1, 'rsa_w']+bb.loc[row-2, 'rsa_w']
                     +bb.loc[row-3, 'rsa_w']+bb.loc[row-4, 'rsa_w']+bb.loc[row-5, 'rsa_w']
                     +bb.loc[row-6, 'rsa_w']+bb.loc[row-7, 'rsa_w']+bb.loc[row-8, 'rsa_w']
                     +bb.loc[row-9, 'rsa_w'])))

J_AVE_p = []
for row in pp.index:
    if pp.loc[row,'seq_w'] != pp.loc[row,'seq_m']:
        J_AVE_p.append(((pp.loc[row, 'delta_rsa']+pp.loc[row+1, 'delta_rsa']+pp.loc[row+2, 'delta_rsa']
                     +pp.loc[row+3, 'delta_rsa']+pp.loc[row+4, 'delta_rsa']+pp.loc[row+5, 'delta_rsa']
                     +pp.loc[row+6, 'delta_rsa']+pp.loc[row+7, 'delta_rsa']+pp.loc[row+8, 'delta_rsa']
                     +pp.loc[row+9, 'delta_rsa']+pp.loc[row-1, 'delta_rsa']+pp.loc[row-2, 'delta_rsa']
                     +pp.loc[row-3, 'delta_rsa']+pp.loc[row-4, 'delta_rsa']+pp.loc[row-5, 'delta_rsa']
                     +pp.loc[row-6, 'delta_rsa']+pp.loc[row-7, 'delta_rsa']+pp.loc[row-8, 'delta_rsa']
                     +pp.loc[row-9, 'delta_rsa'])/19)/((pp.loc[row, 'rsa_w']+pp.loc[row+1, 'rsa_w']+pp.loc[row+2, 'rsa_w']
                     +pp.loc[row+3, 'rsa_w']+pp.loc[row+4, 'rsa_w']+pp.loc[row+5, 'rsa_w']
                     +pp.loc[row+6, 'rsa_w']+pp.loc[row+7, 'rsa_w']+pp.loc[row+8, 'rsa_w']
                     +pp.loc[row+9, 'rsa_w']+pp.loc[row-1, 'rsa_w']+pp.loc[row-2, 'rsa_w']
                     +pp.loc[row-3, 'rsa_w']+pp.loc[row-4, 'rsa_w']+pp.loc[row-5, 'rsa_w']
                     +pp.loc[row-6, 'rsa_w']+pp.loc[row-7, 'rsa_w']+pp.loc[row-8, 'rsa_w']
                     +pp.loc[row-9, 'rsa_w'])/19))  


# In[24]:


sns.distplot(J_AVE_b, hist = True, kde = True, bins = 120, 
                 kde_kws = {'linewidth': 2}, color = 'b')
sns.distplot(J_AVE_p, hist = True, kde = True, bins = 120, 
                 kde_kws = {'linewidth': 2}, color = 'r')

plt.title("Distribution of relative delta rsa values near by variant site")
plt.xlabel("relative rsa values")
plt.ylabel("Frequency")
plt.legend(["benign", "pathogenic"], loc=0)
plt.xlim(-0.2,0.2)


# In[ ]:




