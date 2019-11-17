# -*- coding: utf-8 -*-
"""
Provides source fluxes from the University of `Michigan Radio Observatory 
database <https://dept.astro.lsa.umich.edu/datasets/umrao.php>`_.

The MRO radio source catalog was developed from years of monitoring a select
set of radio sources.  The catalog work was discontinued.

The description and plots included in the UMRAO web site were last updated on 
August 23, 2012. For questions or for additional data, contact 
`Margo Aller <mfa@umich.edu>`_.
"""
import os
import pickle
import urllib

from matplotlib.dates import datestr2num
from numpy import array
from scipy import polyfit, polyval

from Radio_Astronomy import cal_dir
from Data_Reduction import nearest_index

diag = True

pickle_file_path = os.path.join(cal_dir, "michigan_tables.pkl")
pickle_file = open(pickle_file_path,"r")
# pickle_file = open(cal_dir+"michigan_tables.pkl","r")
table_links = pickle.load(pickle_file)
pickle_file.close()

Bnames = table_links.keys()
Bnames.sort()

def get_flux_data(url):
  """
  Gets the latest data from the Michigan database for a source

  The URL is obtained from a local dictionary indexed with the
  source B names.
  """
  table = urllib.urlopen(url)
  table_data = table.read()
  table.close()

  lines = table_data.split('\n')
  num_lines = len(lines)
  times = {}
  fluxes = {}
  sigflux = {}
  for line in lines[3:]:
    if len(line) and line.isspace() == False:
      parts = line.split()
      if len(parts) < 6:
        # No flux data
        print parts
        break
      mjd = float(parts[0])
      date = datestr2num(parts[1])+float(parts[3])/24.
      freq = parts[2]
      flux = float(parts[4])
      sig_flux = float(parts[5])
      if times.has_key(freq) == False:
        times[freq] = []
        fluxes[freq] = []
        sigflux[freq] = []
      times[freq].append(date)
      fluxes[freq].append(flux)
      sigflux[freq].append(sig_flux)
  return times,fluxes,sigflux

def polate_flux(Jname,datenum,freq):
  """
  Interpolate or extrapolate flux of source for given date and frequency
  """
  times, fluxes, sigflux = get_flux_data(table_links[Jname][1])
  # (inter/extra)polate time first
  guess = {}
  for key in times.keys():
    # interpolate or extrapolate?
    datenums = array(times[key])
    if datenum > datenums[-1]:
      # extrapolate forward
      timedelta = datenum - datenums[-1]
      ref_index = nearest_index(datenums,datenums[-1]-timedelta)
      timedata = datenums[ref_index:-1]
      fluxdata = fluxes[key][ref_index:-1]
      (ar,br) = polyfit(timedata,fluxdata,1)
      guess[key] = polyval([ar,br],datenum)
    elif datenum < datenums[0]:
      # extrapolate backward
      time_delta = datenums[0] - datenum
      ref_index = nearest_index(datenums,datenums[0]+timedelta)
      timedata = datenums[0:ref_index]
      fluxdata = fluxes[key][0:ref_index]
      (ar,br)=polyfit(timedata,fluxdata,1)
      guess[key] = polyval([ar,br],datenum)
    else:
      ref_index = nearest_index(datenums,datenum)
      max_index = len(datenums)-1
      first_index = max(ref_index-4,0)
      last_index = min(ref_index+4,max_index)
      timedata = datenums[first_index:last_index]
      fluxdata = fluxes[key][first_index:last_index]
      (ar,br)=polyfit(timedata,fluxdata,1)
      guess[key] = polyval([ar,br],datenum)
  # Now we (inter/extra)polate frequency
  freqs = []
  fluxdata = []
  for key in times.keys():
    freqs.append(float(key))
    fluxdata.append(guess[key])
  (ar,br) = polyfit(freqs,fluxdata,1)
  return polyval([ar,br],freq)

if __name__ == "__main__":
  from matplotlib.dates import datestr2num
  from pylab import *
  from Radio_Astronomy.vla_cal import get_cal_data, get_cal_dict, \
                                      Jnames_to_B, VLA_name_xref
  cal_data = get_cal_dict()
  vla_sources = cal_data.keys()
  Bnames_dict,cat_3C_dict = VLA_name_xref(cal_data)
  Jnames_dict = Jnames_to_B(Bnames_dict)
  freqs = [4.8, 8.0, 14.5]
  source = "3C84"
  Jname = cat_3C_dict[source]
  Bname = Jnames_dict[Jname]
  flux = polate_flux(Bname,datestr2num("2012 Jul 4 00:00 UT"),freqs)
  plot(freqs,flux,'.')
  x = linspace(1.,30.)
  y = polate_flux(Bname,datestr2num("2012 Jul 4 00:00 UT"),x)
  plot(x,y,'-')
  grid()
  show()
