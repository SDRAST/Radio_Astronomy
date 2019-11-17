"""
conversion from waveguide band to frequency and vice-versa

Examples::

 In [1]: from Radio_Astronomy.bands import *
 In [2]: frequency_to_band(3)
 Out[2]: 'S'
 In [3]: band_to_frequency('K')
 Out[3]: 22
"""
def frequency_to_band(freq):
  """
  band code from frequency in GHz
  """
  if                  freq  <  1:   return None
  elif freq >= 1   and freq <  2:   return "L"
  elif freq >= 2   and freq <  4:   return "S"
  elif freq >= 4   and freq <  8:   return "C"
  elif freq >= 8   and freq < 12:   return "X"
  elif freq >=12   and freq < 18:   return "Ka"
  elif freq >=18   and freq < 26.5: return "K"
  elif freq >=26.5 and freq < 40:   return "Ka"
  elif freq >=40   and freq < 50:   return "Q"
  elif freq >=50   and freq < 75:   return "V"
  elif freq >=75   and freq <115:   return "W"
  else:                             return "D"

def band_to_frequency(band):
  """
  nominal band center frequency in GHz from band code
  """
  if   band == "L":  return 1.7
  elif band == "S":  return 2.3
  elif band == "C":  return 5.0
  elif band == "X":  return 8.45
  elif band == "Ku": return 15
  elif band == "K":  return 22
  elif band == "Ka": return 34
  elif band == "Q":  return 42
  elif band == "W":  return 90
  else: return None
