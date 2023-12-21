import requests

def get_file_from(url):
  """Retrieves file from http address
  :param str address: the http address of the file
  :return str file: the content of the file
  :return None: if the request fails, returns None"""
  r = requests.get(url)
  if r.status_code == 200:
    return r.text
  return None