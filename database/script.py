import os
import api
import parse_summary as ps
import re

# online address of NCBI database
archea_address = "http://ftp.ncbi.nih.gov/genomes/refseq/archaea/"
bacteria_address = "http://ftp.ncbi.nih.gov/genomes/refseq/bacteria/"

# local path to datas
archea_path = "./database/archea/"
bacteria_path = "./database/bacteria/"
if not os.path.exists(archea_path):
  os.mkdir(archea_path)
if not os.path.exists(bacteria_path):
  os.mkdir(bacteria_path)

# 1. Retrieve data summaries
summary_file = "assembly_summary.txt"
if not os.path.exists(bacteria_path + summary_file):
  bacteria_summary = api.get_file(bacteria_address + summary_file, bacteria_path + summary_file)

if not os.path.exists(archea_path + summary_file):
  archea_summary = api.get_file(archea_address + summary_file, archea_path + summary_file)

# 2. Get list of whole genomes
bacteria_data = ps.read_summary(bacteria_path + summary_file)
bacteria_dict = ps.select_full_genome(bacteria_data)
archea_data = ps.read_summary(archea_path + summary_file)
archea_dict = ps.select_full_genome(archea_data)

# 3. Retrieve .fna of each species in list
for species in archea_dict.keys():
  file = archea_path + re.sub('[^0-9a-zA-Z]+', '_', species) + '.fna.gz'
  if not os.path.exists(file):
    print(f'Downloading {species}')
    api.get_file(archea_dict[species], file)
  
for species in bacteria_dict.keys():
  file = bacteria_path + re.sub('[^0-9a-zA-Z]+', '_', species) + '.fna.gz'
  if not os.path.exists(file):
    print(f'Downloading {species}')
    api.get_file(bacteria_dict[species], file)