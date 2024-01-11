import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def stream_kmers(file, k):
  """ Transforms a file of nucleotid sequences into its binary kmer representation
  :param Array file: file of sequences
  :param int k: the number of nucleotides involved into each kmer
  :return generator kmer: the integer representation of each of the kmers
  """
  for seq in file:
    kmer = 0
    kmer_rev_comp = 0
    mask = (1<<(2*(k-1))) - 1
    for i in range(len(seq)):
      nuc = ((ord(seq[i]) >> 1) & 0b11)
      kmer <<= 2
      kmer = kmer + nuc
      if i >= (k - 1):
        # computing complementary nucleotid 
        nuc_comp = ((nuc + 0b10) & 0b11) << ((k - 1)*2)
        kmer_rev_comp += nuc_comp
        # we return the minimum of reverse complement and original sequence to keep the smallest key
        yield min(kmer, kmer_rev_comp)
        # we delete the data
        kmer_rev_comp >>= 2
        kmer &= mask
      else:
        # computing complementary nucleotid 
        nuc_comp = ((nuc + 0b10) & 0b11) << (i*2)
        kmer_rev_comp += nuc_comp

def kmer2str(val, k):
  """ Transform a kmer integer into a its string representation
  :param int val: An integer representation of a kmer
  :param int k: The number of nucleotides involved into the kmer.
  :return str: The kmer string formatted
  """
  letters = ['A', 'C', 'T', 'G']
  str_val = []
  for i in range(k):
    str_val.append(letters[val & 0b11])
    val >>= 2

  str_val.reverse()
  return "".join(str_val)

def count_kmers(sequences, k):
  """Counts the occurences of each kmer in sequences
  :param list sequences: a list of sequences
  :param int k: kmer length
  :return dict res: dictionary of binary kmers occurences
  """
  res = {}
  for kmer in stream_kmers(sequences, k):
    if kmer not in res:
      res[kmer] = 1
    else:
      res[kmer] += 1
  return res

def count_aa_kmers(sequences, k):
  """Counts the occurences of each kmer in proteic sequences
  :param list sequences: a list of sequences
  :param int k: kmer length
  :return dict res: dictionary of binary kmers occurences
  """
  res = {}
  for kmer in stream_aa_kmers(sequences, k):
    if kmer not in res:
      res[kmer] = 1
    else:
      res[kmer] += 1
  return res

def plot_spectrum(counts, name, k, saveto, bins=100):
  """Plots kmer spectrum and saves in saveto file
  :param Array counts: occurrences of the kmers
  :param str name: species name
  :param int k: kmer length
  :param str saveto: path to the output file
  :param int bins: number of bins to consider in the histogram"""
  plt.hist(counts, range = (0, max(counts)), bins = bins)
  plt.title(f'{name} {k}-mer spectrum')
  plt.xlabel('Occurence')
  plt.ylabel('Frequence')
  plt.savefig(saveto)
  plt.close()

def write_profile(counts, outfile, k):
  """Writes counts as csv in outfile
  :param dict counts: dicionnary of occurence per kmer (binary)
  :param str outfile: path to outfile
  :param int k: length of kmer"""
  with open(outfile, "w") as f:
    f.write('Kmer\tBinary kmer\tCount\n')
    for kmer in counts.keys():
      f.write(f'{kmer2str(kmer, k)}\t{kmer}\t{counts[kmer]}\n')


def profile_fusion2(profile_list, path, k, outfile):
  """Creates a table fusioning all profiles
  :param Array profile_list: list of files to fuse
  """
  total_counts = {}
  species_list = []

  # reading all profiles
  for file in profile_list:
    print(file)
    species = file[:-len('_profile.txt')]
    species_list += [species]
    data = pd.read_csv(path + file, sep='\t')
    for _, row in data.iterrows():
      kmer = row['Binary kmer']
      if kmer not in total_counts:
        total_counts[kmer] = [(species,row['Count'])]
      else:
        total_counts[kmer] += [(species,row['Count'])]

  bin_kmer_list = list(total_counts.keys())
  nb_species = len(species_list)

  # writting profiles dictionnary to table
  table = []
  for kmer in bin_kmer_list:
    row = [0] * nb_species
    for sp, count in total_counts[kmer]:
      index = species_list.index(sp)
      row[index] = count
    table += [row]

  table = pd.DataFrame(np.transpose(table), index = species_list, columns=bin_kmer_list)
  table.to_csv(outfile, sep='\t')

def profile_fusion(profile_list, path, k, outfile):
  """Creates a table fusioning all profiles
  :param Array profile_list: list of files to fuse
  """
  data = pd.DataFrame([])

  # reading all profiles
  for file in profile_list:
    print(file)
    data_file = pd.read_csv(path + file, sep='\t')
    species = file[:-len('_profile.txt')]
    values = pd.DataFrame(np.array(data_file['Count']).reshape(1, -1), index=[species], columns=data_file['Binary kmer'])
    data = pd.concat([data, values])

  data.fillna(0)
  data.to_csv(outfile, sep='\t')

# TEST: the input sequence, put through the two methods, must be found in output
seqs=['AGTCTGAGTGC']
k = len(seqs[0])
values = stream_kmers(seqs, k)
out_seq=''
for val in values:
  out_seq = kmer2str(val, k)
assert(seqs[0] == out_seq)

seqs=['AAAA']
k = len(seqs[0])
values = stream_kmers(seqs, k)
outseq=''
for val in values:
  out_seq = kmer2str(val, k)
assert(seqs[0] == out_seq)


dico_aa = {
    "A":0,
    "F":1,
    "V":2,
    "G":3,
    "L":4,
    "I":5,
    "T":6,
    "S":7,
    "Y":8,
    "N":9,
    "D":10,
    "E":11,
    "H":12,
    "K":13,
    "R":14,
    "M":15,
    "W":16,
    "Q":17,
    "C":18,
    "P":19
}

def kmer2str_aa(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ["A","F","V","G","L","I","T","S","Y","N","D","E","H","K","R","M","W","Q","C","P"]
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11111])
        val >>= 5
    str_val.reverse()
    return "".join(str_val)

def stream_aa_kmers(sequences, k):
    x = 0
    mask = (1<<(k*5))-1
    for seq in sequences:
      for aa in range(0,len(seq)):
          if seq[aa] not in dico_aa:
              pass
          else:
              n = dico_aa[seq[aa]]&0b11111
              x <<=5
              x += n
              x&=mask
              if aa > k-2:
                  yield x
