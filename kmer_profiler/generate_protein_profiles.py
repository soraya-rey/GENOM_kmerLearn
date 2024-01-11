from read_files import read_fna
from kmer import *
import os
from multiprocessing import Process
import random

archea_database = './database/archea_protein/'
bacteria_database = './database/bacteria_protein/'

archea_profiles = './profiles/archea_protein/'
bacteria_profiles = './profiles/bacteria_protein/'

# 1. Create profiles destination folders
if not os.path.exists('./profiles/'):
  os.mkdir('./profiles/')

if not os.path.exists(archea_profiles):
  os.mkdir(archea_profiles)
if not os.path.exists(bacteria_profiles):
  os.mkdir(bacteria_profiles)

# 2. Retrieve database files 
archea_files = [file for file in os.listdir(archea_database) if file[-7:]== '.fna.gz']
bacteria_files = [file for file in os.listdir(bacteria_database) if file[-7:]== '.fna.gz']

# 3. For each whole genome fasta, create kmer profile for several values of k
k_values = {3: 25, 5: 100, 8:100, 6:100}

for k in k_values.keys():
    path = archea_profiles + f'{k}-mers_profiles/'
    if not os.path.exists(path):
      os.mkdir(path)

    path = bacteria_profiles + f'{k}-mers_profiles/'
    if not os.path.exists(path):
      os.mkdir(path)

# For archeas
nb = len(archea_files)
i=0
for file in archea_files:
  i+=1
  name = file.split('/')[-1].split('.')[0]
  # checking if profiles already exists
  k_to_do = [k for k in k_values.keys() if not os.path.exists(archea_profiles + f'{k}-mers_genes_profiles/' + f'{name}_profile.txt')]
  if len(k_to_do) == 0:
    continue

  # reading sequences in fasta
  dict = read_fna(archea_database + file)
  print(f'Computing profiles for {name} ({round((i/nb)*100, 2)}% : {i}/{nb})')

  # for each value of k, compute k-mer profile
  for k in k_to_do:
    path = archea_profiles + f'{k}-mers_profiles/'
    outfile_path = path + name

    for index in np.array(list(dict.keys()))[random.sample(range(len(dict)), 50)]:
      # counting kmers
      counts = count_aa_kmers([dict[index]], k)
      occurences = counts.values()
      protein = index.split('.')[0]

      # writting profile spectrum and counts
      write_profile(counts, f'{outfile_path}.{protein}_profile.txt', k)

# For bacterias
nb = len(bacteria_files)
i = 0
for file in bacteria_files:
  i += 1
  name = file.split('/')[-1].split('.')[0]

  # checking if profiles already exists
  k_to_do = [k for k in k_values.keys() if not os.path.exists(bacteria_profiles + f'{k}-mers_profiles/' + f'{name}_profile.txt')]
  if len(k_to_do) == 0:
    continue

  # reading sequences in fasta
  dict = read_fna(bacteria_database + file)
  print(f'Computing profiles for {name} ({round((i/nb)*100, 2)}% : {i}/{nb})')

  # for each value of k, compute k-mer profile
  for k in k_values.keys():
    path = bacteria_profiles + f'{k}-mers_profiles/'

    # reading sequences in fasta
    outfile_path = path + name

    # counting kmers
    for index in np.array(list(dict.keys()))[random.sample(range(len(dict)), 50)]:
      counts = count_aa_kmers([dict[index]], k)
      occurences = counts.values()
      protein = index.split('.')[0]

      # writting profile spectrum and counts
      write_profile(counts, f'{outfile_path}.{protein}_profile.txt', k)


# computing unique files combining all profiles for k-mers 
for k in k_values.keys():
  outfile = f'{k}-mer_archea_protein.csv'
  path = archea_profiles + f'{k}-mers_profiles/'
  if not os.path.exists('./profiles/' + outfile):
    files = [file for file in os.listdir(path) if file[-len('.txt'):]=='.txt']
    profile_fusion(files, path, k, './profiles/' + outfile)

  outfile = f'{k}-mer_bacteria_protein.csv'
  path = bacteria_profiles + f'{k}-mers_profiles/'
  if not os.path.exists('./profiles/' + outfile):
    files = [file for file in os.listdir(path) if file[-len('.txt'):]=='.txt']
    profile_fusion(files, path, k, './profiles/' + outfile)
