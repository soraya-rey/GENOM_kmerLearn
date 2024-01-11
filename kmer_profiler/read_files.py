import gzip
from Bio import SeqIO

def read_fna(filepath):
  """Reads a fasta file of a complete genome sequence
  :param str filepath: path to the fasta file (format: .fna.gz)
  :return str name: species name
  :return dict sequences: dictionnary of the sequences indexed by ids"""
  sequences = {}
  with gzip.open(filepath, 'rt') as handle:
    for record in SeqIO.parse(handle, "fasta"):
      sequences[record.id] = str(record.seq)
  return sequences