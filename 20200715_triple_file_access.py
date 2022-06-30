#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import pickle as pk
from matplotlib import pyplot as plt


# BEFORE RUNNING, replace `[PATH]` in the following cell AND in `./C3POa/config.txt` with the full path to the `cas9-random-access` cloned repo.

# In[2]:


date = "20200715"
input_dir = "/[PATH]/cas9-random-access/%s_basecalled_reads" % date # absolute path of directory containing basecalled reads in fastq format
output_dir = "/[PATH]/cas9-random-access" # absolute path of output directory
nuc_mat = "/[PATH]/cas9-random-access/C3POa/NUC.4.4.mat"
splint_seq = "splint_all.fasta"
config = "C3POa/config.txt"


# In[3]:


reads_basecalled = 0
for x in os.listdir(input_dir):
    fastq = open(os.path.join(input_dir, x), 'r')
    for line in fastq:
        if line.strip() == "+":
            reads_basecalled += 1
    fastq.close()
print("Total Reads Passed Basecalling:", reads_basecalled)


# # Preprocessing
# Demultiplex reads based on their file addresses.

# In[3]:


pre_dir = os.path.join(output_dir, "preprocessing") # store demultiplexed reads in this directory


# In[4]:


q_score = 9 
pre_length_cutoff = 100 


# In[6]:


for fastq in os.listdir(input_dir):
    print(fastq)
    fastq = os.path.join(input_dir, fastq)
    preprocess_cmd = "C3POa/C3POa_preprocessing.py -i %s -o %s -q %d -l %d -s %s -c %s" % (fastq, pre_dir, q_score, pre_length_cutoff, splint_seq, config)
    get_ipython().system('python3 {preprocess_cmd}')
    print()


# In[7]:


reads_per_guide = dict.fromkeys(range(1,27), 0)
for x in os.listdir(pre_dir):
    for guide in range(1,27):
        if "splint_g%02d" % guide in os.listdir(os.path.join(pre_dir, x)):
            splint_dir = os.path.join(pre_dir, x, "splint_g%02d" % guide)
            if os.listdir(splint_dir):
                fastq = open(os.path.join(splint_dir, "R2C2_raw_reads.fastq"), 'r')
                for line in fastq:
                    if line.strip() == "+":
                        reads_per_guide[guide] += 1


# In[8]:


print("Total Reads After Preprocessing:", sum(reads_per_guide.values()))
print("% of Basecalled Reads:", 100. * sum(reads_per_guide.values()) / reads_basecalled)


# # Processing
# Consolidate repeats in each concatemer into one consensus read using C3POa. 

# In[5]:


consensus_dir = os.path.join(output_dir, "consensus") # store consensus reads in this directory
temp_dir = os.path.join(consensus_dir, "temp")


# In[6]:


raw_seq_length_cutoff = 100
peak_dist_cutoff = 100 # Median distance between peaks cutoff. This should be the length of your shortest input sequence in your library preparation.


# In[ ]:


for d in os.listdir(pre_dir):
    for guide in range(1,27):
        reads_file = os.path.join(pre_dir, d, "splint_g%02d" % guide, "R2C2_raw_reads.fastq")
        if os.path.isfile(reads_file):
            current_temp_dir = os.path.join(temp_dir, "temp_dir%s_g%02d" % (d, guide))
            if not os.path.exists(current_temp_dir):
                os.makedirs(current_temp_dir)
            out_file = os.path.join(consensus_dir, "R2C2_consensus_g%02d.fasta" % (guide))
            partial_reads_file = os.path.join(consensus_dir, "R2C2_partial_reads_g%02d.fasta" % (guide))
            process_cmd = "C3POa/C3POa.py -r %s -p %s -m %s -l %d -d %d -c %s -o %s -s %s" % (reads_file, current_temp_dir, nuc_mat, raw_seq_length_cutoff, peak_dist_cutoff, config, out_file, partial_reads_file)
            get_ipython().system('python3 {process_cmd}')


# In[8]:


copies_per_guide = dict.fromkeys(range(1,27))
for x in copies_per_guide:
    copies_per_guide[x] = []
for guide in range(1,27):
    fasta = open(os.path.join(consensus_dir, "R2C2_consensus_g%02d.fasta" % guide), 'r')
    for line in fasta:
        if line[0] == ">":
            copies_per_guide[guide].append(int(line.split("_")[3]) + 1)
    fasta.close()
    fasta = open(os.path.join(consensus_dir, "R2C2_partial_reads_g%02d.fasta" % guide), 'r')
    for line in fasta:
        if line[0] == ">":
            copies_per_guide[guide].append(1)
    fasta.close()


# In[11]:


print("Total Reads After Processing:", sum([len(copies_per_guide[x]) for x in copies_per_guide]))
print("% of Basecalled Reads:", 100. * sum([len(copies_per_guide[x]) for x in copies_per_guide]) / reads_basecalled)


# # Analysis

# In[12]:


normalization_dict = pk.load(open("normalization_dict.pkl", "rb"))


# In[13]:


fig, ax = plt.subplots(figsize=(8,4))
for guide, x in copies_per_guide.items():
    if guide == 26:
        continue
    if guide in [2,13,24]:
        ax.bar(guide, float(len(x)) / sum([len(copies_per_guide[x]) for x in copies_per_guide]), color="#6EB3E4")
    else:
        ax.bar(guide, float(len(x)) / sum([len(copies_per_guide[x]) for x in copies_per_guide]), color="silver")
        
plt.xticks(list(copies_per_guide.keys())[:-1])
ax.set_yscale('log')
plt.xlabel("File #")
plt.ylabel("Fraction of Total Concatemers")
plt.title("%s - Concatemers Per File" % date)
plt.show()


# In[14]:


fig, ax = plt.subplots(figsize=(6,4))
for guide in range(1,26):
    if guide in [2,13,24]:
        ax.scatter([1], float(len(copies_per_guide[guide])) / sum([len(x) for x in copies_per_guide.values()]), color='purple', marker='x')
    else:
        ax.scatter([2], float(len(copies_per_guide[guide])) / sum([len(x) for x in copies_per_guide.values()]), color='purple', marker='x')
ax.set_xticks(range(1,3))
ax.set_xticklabels(['Accessed', 'Unaccessed'])
ax.set_yscale('log')
plt.xlim(0.5,2.5)
plt.ylabel("Fraction of Total Concatemers")
plt.title('%s - Fraction of Concatemers Read for Accessed vs. Unaccessed Files' % date)
plt.show()


# In[15]:


fig, ax = plt.subplots(figsize=(8,4))
for guide, x in copies_per_guide.items():
    if guide == 26:
        continue
    if guide in [2,13,24]:
        plt.bar(guide, float(len(x)) / sum([len(copies_per_guide[x]) for x in copies_per_guide]) / normalization_dict[guide],
                color="#6EB3E4")
    else:
        plt.bar(guide, float(len(x)) / sum([len(copies_per_guide[x]) for x in copies_per_guide]) / normalization_dict[guide], 
                color="silver")
        
plt.xticks(list(copies_per_guide.keys())[:-1])
ax.set_yscale('log')
plt.xlabel("File #")
plt.ylabel("Enrichment Score")
plt.title("%s - Concatemers Per File (Enrichment Scores)" % date)
plt.show()


# In[16]:


fig, ax = plt.subplots(figsize=(6,4))
for guide in range(1,26):
    if guide in [2,13,24]:
        ax.scatter([1], float(len(copies_per_guide[guide])) / sum([len(x) for x in copies_per_guide.values()]) / normalization_dict[guide], color='purple', marker='x')
    else:
        ax.scatter([2], float(len(copies_per_guide[guide])) / sum([len(x) for x in copies_per_guide.values()]) / normalization_dict[guide], color='purple', marker='x')
ax.set_xticks(range(1,3))
ax.set_xticklabels(['Accessed', 'Unaccessed'])
ax.set_yscale('log')
plt.xlim(0.5,2.5)
plt.ylabel("Enrichment Score")
plt.title('%s - Enrichment Score for Accessed vs. Unaccessed Files' % date)
plt.show()


# In[18]:


plt.figure(figsize=(8,4))
for guide, x in copies_per_guide.items():
    if guide == 26:
        continue
    if guide in [2,13,24]:
        plt.bar(guide, np.mean(x), color="#6EB3E4")
    else:
        plt.bar(guide, np.mean(x), color="silver")
plt.xticks(list(copies_per_guide.keys())[:-1])
plt.xlabel("File #")
plt.ylabel("Mean # of Copies per Concatemer")
plt.title("%s - Number of Copies Per Concatemer" % date)
plt.show()


# In[19]:


fig, ax = plt.subplots(figsize=(6,4))
for guide in range(1,26):
    if guide in [2,13,24]:
        ax.scatter([1], np.mean(copies_per_guide[guide]), color='purple', marker='x')
    else:
        ax.scatter([2], np.mean(copies_per_guide[guide]), color='purple', marker='x')
ax.set_xticks(range(1,3))
ax.set_xticklabels(['Accessed', 'Unaccessed'])
plt.xlim(0.5,2.5)
plt.ylabel("# of Copies per Concatemer")
plt.title('%s - Mean Concatemer Length for Accessed vs. Unaccessed Files' % date)
plt.show()

