import os
import api
import parse_summary as ps

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
  bacteria_summary = api.get_file_from(bacteria_address + summary_file)
  with open(bacteria_path + summary_file, 'w', encoding="utf-8") as file:
    file.write(bacteria_summary)
if not os.path.exists(archea_path + summary_file):
  archea_summary = api.get_file_from(archea_address + summary_file)
  with open(archea_path + summary_file, 'w', encoding="utf-8") as file:
    file.write(archea_summary)

# 2. Get list of whole genomes
bacteria_data = ps.read_summary(bacteria_path + summary_file)
bacteria_dict = ps.select_full_genome(bacteria_data)
archea_data = ps.read_summary(archea_path + summary_file)
archea_dict = ps.select_full_genome(archea_data)

# 3. Retrieve .fna of each species in list
for species in archea_dict.keys():
  file = archea_path + species.replace(' ', '_').replace('/', '-') + '.fna.gz'
  if not os.path.exists(file):
    data = api.get_file_from(archea_dict[species])
    with open(file, 'w', encoding="utf-8") as f:
      f.write(data)
  
for species in bacteria_dict.keys():
  file = bacteria_path + species.replace(' ', '_').replace('/', '-') + '.fna.gz'
  if not os.path.exists(file):
    data = api.get_file_from(bacteria_dict[species])
    with open(file, 'w', encoding="utf-8") as f:
      f.write(data)