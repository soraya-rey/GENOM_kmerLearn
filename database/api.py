from urllib.request import urlretrieve

def get_file(url, filepath):
  """Retrieves file from url address and saves it un filepath
  :param str address: the http address of the file
  :return str file: the content of the file
  :return None: if the request fails, returns None"""
  urlretrieve(url, filepath)