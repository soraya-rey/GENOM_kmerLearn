import numpy as np
import pandas as pd

def read_summary(file):
  """Reads a summary file
  :param str file: path to the summary file
  :return DataFrame data: content of the file as an pandas DataFrame array"""
  table = []
  columns = []
  with open(file,'r', encoding="utf-8") as f:
    for line in f.readlines():
      # line starting with # is columns line
      if line[0] =="#":
        if line[1] != "#":
          columns = line.replace("\n", "").replace("#", "").split("\t")
      else:
        table += [line.replace("\n", "").split("\t")]
  data = pd.DataFrame(table, columns = columns)
  return data

def select_full_genome(data, list):
  """Returns names of species for which the genome is indicated as full in a summary file and is present in list
  :param DataFrame data: content of the file as an pandas DataFrame array
  :return dictionary genome_list: dictionary of species with full genome (name:address)"""
  full_lines = data[data["assembly_level"] == "Complete Genome"]
  genome_list = {}
  for _, line in full_lines.iterrows():
    organism = line["organism_name"]
    if organism in list:
      path = line["ftp_path"] + "/" + line["ftp_path"].split("/")[-1] + "_genomic.fna.gz"
      genome_list[organism] = path
  return genome_list

def select_cds(data, list):
  """Returns names of species for which the genome is indicated as full in a summary file
  :param DataFrame data: content of the file as an pandas DataFrame array
  :return dictionary genome_list: dictionary of species with full genome (name:address)"""
  full_lines = data[data["assembly_level"] == "Complete Genome"]
  genome_list = {}
  for _, line in full_lines.iterrows():
    organism = line["organism_name"]
    if organism in list:
      path = line["ftp_path"] + "/" + line["ftp_path"].split("/")[-1] + "_cds_from_genomic.fna.gz"
      genome_list[organism] = path
  return genome_list

def select_prot(data, list):
  """Returns names of species for which the genome is indicated as full in a summary file
  :param DataFrame data: content of the file as an pandas DataFrame array
  :return dictionary genome_list: dictionary of species with full genome (name:address)"""
  full_lines = data[data["assembly_level"] == "Complete Genome"]
  genome_list = {}
  for _, line in full_lines.iterrows():
    organism = line["organism_name"]
    if organism in list:
      path = line["ftp_path"] + "/" + line["ftp_path"].split("/")[-1] + "_protein.faa.gz"
      genome_list[organism] = path
  return genome_list

def read_list(file):
  """Reads a list of species in txt file"""
  results = []
  with open(file) as f:
    lines = f.readlines()
    for line in lines:
      results += [line.replace('\n', '')]
  return results