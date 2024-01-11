import os
import api
import parse_summary as ps
import re

# online address of NCBI database
archea_address = "http://ftp.ncbi.nih.gov/genomes/refseq/archaea/"
bacteria_address = "http://ftp.ncbi.nih.gov/genomes/refseq/bacteria/"

# local path to datas
archea_path = "./database/archea/"
archea_cds_path = "./database/archea_cds/"
archea_protein_path = "./database/archea_protein/"
bacteria_path = "./database/bacteria/"
bacteria_cds_path = "./database/bacteria_cds/"
bacteria_protein_path = "./database/bacteria_protein/"

archea_list_path = './database/archea_selection.txt'
bacteria_list_path = './database/bacteria_selection.txt'
archea_list = ps.read_list(archea_list_path)
bacteria_list = ps.read_list(bacteria_list_path)

if not os.path.exists(archea_path):
  os.mkdir(archea_path)
if not os.path.exists(bacteria_path):
  os.mkdir(bacteria_path)
if not os.path.exists(archea_cds_path):
  os.mkdir(archea_cds_path)
if not os.path.exists(bacteria_cds_path):
  os.mkdir(bacteria_cds_path)
if not os.path.exists(archea_protein_path):
  os.mkdir(archea_protein_path)
if not os.path.exists(bacteria_protein_path):
  os.mkdir(bacteria_protein_path)

# 1. Retrieve data summaries
summary_file = "assembly_summary.txt"
if not os.path.exists(bacteria_path + summary_file):
  bacteria_summary = api.get_file(bacteria_address + summary_file, bacteria_path + summary_file)

if not os.path.exists(archea_path + summary_file):
  archea_summary = api.get_file(archea_address + summary_file, archea_path + summary_file)

# 2. Get list of whole genomes, cds and protein files to download
bacteria_data = ps.read_summary(bacteria_path + summary_file)
bacteria_dict = ps.select_full_genome(bacteria_data, bacteria_list)
bacteria_cds_dict = ps.select_cds(bacteria_data, bacteria_list)
bacteria_protein_dict = ps.select_prot(bacteria_data, bacteria_list)
archea_data = ps.read_summary(archea_path + summary_file)
archea_dict = ps.select_full_genome(archea_data, archea_list)
archea_cds_dict = ps.select_cds(archea_data, archea_list)
archea_protein_dict = ps.select_prot(archea_data, archea_list)

# 3. Retrieve genome.fna of each species in list
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

# 4. Retrieve cds of each species in list
for species in archea_cds_dict.keys():
  file = archea_cds_path + re.sub('[^0-9a-zA-Z]+', '_', species) + '.fna.gz'
  if not os.path.exists(file):
    print(f'Downloading {species}')
    api.get_file(archea_cds_dict[species], file)
  
for species in bacteria_cds_dict.keys():
  file = bacteria_cds_path + re.sub('[^0-9a-zA-Z]+', '_', species) + '.fna.gz'
  if not os.path.exists(file):
    print(f'Downloading {species}')
    api.get_file(bacteria_cds_dict[species], file)

# 5. Retrieve proteins .faa of each species in list
for species in archea_protein_dict.keys():
  file = archea_protein_path + re.sub('[^0-9a-zA-Z]+', '_', species) + '.fna.gz'
  if not os.path.exists(file):
    print(f'Downloading {species}')
    api.get_file(archea_protein_dict[species], file)
  
for species in bacteria_protein_dict.keys():
  file = bacteria_protein_path + re.sub('[^0-9a-zA-Z]+', '_', species) + '.fna.gz'
  if not os.path.exists(file):
    print(f'Downloading {species}')
    api.get_file(bacteria_protein_dict[species], file)
