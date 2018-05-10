# -*- coding: utf-8 -*-
"""
module vla_cal

Module vla_cal creates or updates files of VLA calibrators and retrieves
data for a specified calibrator.

Administration
==============
It is expected that the pickle files VLA_cals, B_names and 3C_names exist.
If they do not, the following must be done in the module directory
/usr/local/lib/python2.6/site-packages/Radio_Astronomy/:
>>> from Radio_Astronomy.vla_cal import *
>>> cal_data = get_VLA_calibrators()
>>> make_VLA_pickle_file(cal_data)
>>> Bname_dict,cat_3C_dict = VLA_name_xref(cal_data)
>>> make_Bname_pickle_file(Bname_dict)
>>> make_3C_pickle_file(cat_3C_dict)

3C_VLA_cals is a convenience.pickle file that provides coordinate data
keyed to 3C names.  It is created as follows:
>>> make_3C_file()

User Instructions
=================
If you only needs coordinates for a 3C source then this will suffice:
>>> ra,dec,jname = get_3C_coords(name)

If you need all the data for a source then do this:
>>> data_dict = get_cal_data(source)

Notes
=====
One source of annoyance is that astronomers do not all or always use the same
designations as NRAO, that is, more or less characters, like b1937+21,
j0218+4232 and so on.

IAU source designation recommendations can be found at
http://cdsweb.u-strasbg.fr/Dic/iau-spec.html

"""
diag = False

import logging
import urllib
import pickle
import re
from math import pi
from Astronomy import formats
import struct
from scipy import polyfit, polyval

from Radio_Astronomy import cal_dir

module_logger = logging.getLogger(__name__)

def make_3C_file(url='http://www.vla.nrao.edu/astro/calib/manual/csource.html'):
    """
    Makes VLA calibrators 3C pickle file.

    The 3C names are the keys.  The data for each key is a list with formatted
    J2000 right ascension and declination, and the IAU name.

    @param url : location of the VLA calibrator database
    @type  url : str

    @return: True
    """
    vla_cal = \
            urllib.urlopen(url)
    file_open = True
    cal_data = {}
    while file_open:
        line = vla_cal.readline()
        # get the lines with J2000 coordinates and the 3C name
        if re.search('J2000',line[0:20]):
            source_data = line.split()
            if len(source_data) == 7:
                if re.match('3C',source_data[6]):
                    jname,epoch,type,ra,dec,date,name = source_data
                    data = [ra,dec,jname]
                    cal_data[name] = data
        if line == '':
            vla_cal.close()
            file_open = False
    dbfile = open(cal_dir+'3C_VLA_cals','w')
    pickle.dump(cal_data,dbfile)
    dbfile.close()
    return True

def get_3C_coords(name):
    """
    Formatted J2000 right ascension and declination and IAU name

    Returns the formatted J2000 right ascension and declination and IAU name
    given the 3C name.
    Example:
    >>> ra,dec,iau = get_3C_coords('3C286')
    >>> print ra,dec,iau
    13h31m08.287984s 30d30'32.958850" 1331+305

    @param name : 3C name, like 3C123

    @return: ra, dec, IAU_name
    """
    dbfile = open(cal_dir+'3C_VLA_cals','r')
    data = pickle.load(dbfile)
    dbfile.close()
    return data[name]

def get_VLA_calibrators(url='http://www.vla.nrao.edu/astro/calib/manual/csource.html'):
  """
  Create a dictionary keyed on J-names of calibrators in the VLA data base.

  Entries in the VLA database are separated by spaces. This is a typical entry::
   0001+192   J2000  A 00h01m08.621563s  19d14'33.801860"  Aug01  JVAS
   2358+189   B1950  A 23h58m34.865400s  18d57'51.753000"
   -----------------------------------------------------
   BAND        A B C D    FLUX(Jy)    UVMIN(kL)  UVMAX(kL)
   =====================================================
   0.7cm    Q  W W W W       0.18

   These data are converted to this format:
   cals['0001+192'] = {'bname': '2358+189', 'dec': 19.24272273888889,
                       'mm7': 0.17999999999999999, 'ra': 0.019061545277777779}

  @param url : URL of the VLA calibrator source list

  @return: dictionary of dictionaries
  """

  f = urllib.urlopen(url)
  # Skip over title, internal links, and column headings
  for i in range(11):
    line = f.readline()
  # Pull in the data
  catalog = f.readlines()
  f.close()
  n_lines = len(catalog)

  cal_data = {}
  get_fluxes = False # Get source name and coordinates first
  for i in range(n_lines):
    line = catalog[i].strip('\n')
    # This strips out any HTML tags
    line,count = re.subn('<.*?>','',line)
    if line == '' and count == 1:
      # The last line is '</PRE>' so it becomes an empty line and we are done
      cal_data[current_source] = this_source
      break
    elif line.isspace():
      # A blank line starts reading a new source block.
      # Start a new dictionary for this source
      try:
        # If there is a dictionary with data for a source, put it in the
        # cal_data dictionary for this source and start a new source dictionary
        if len(this_source) > 0:
          cal_data[current_source] = this_source
      except NameError:
        # No source has yet been processed
        pass
      this_source = {}
      # We are not yet looking for flux data
      get_fluxes == False
    elif re.search('J2000',line):
      # This is the first line of a source block.  It has the J-name which is
      # a key in cal_data
      data = line.split()
      current_source = data[0]
      cal_data[current_source] = {}
      try:
        this_source['ra'] = \
                        formats.hms_delimited_angle_to_rads(data[3])*12/pi
      except Exception, details:
        # This should not happen; quit and diagnose manually
        module_logger.error("Could not parse RA:", data[3])
        module_logger.error('in:',line)
        module_logger.error(details)
        break
      try:
        this_source['dec'] = \
                        formats.dms_delimited_angle_to_rads(data[4])*180/pi
      except Exception, details:
        # This should not happen; quit and diagnose manually
        module_logger.error("Could not parse decl.:", data[4])
        module_logger.error("in:",line)
        module_logger.error(details)
        break
      # check for alternate names
      altname = line[63:].strip()
      if re.search('3C',altname):
        this_source['cat3c']=altname.strip()
      else:
        if not altname.isspace():
          module_logger.debug(current_source,"has alternate name",altname,", length",len(altname))
    elif re.search('B1950',line):
      # This is the B1950 name; don't bother to get the 1950 coordinates
      data = line.split()
      this_source['bname'] = data[0]
    elif line[:2] == '--':
      # This indicates the end of the coordinates subsection
      pass
    elif line[:4] == 'BAND':
      # This is the header of the flux subsection
      pass
    elif line[:2] == '==':
      # This is the separator between header and data of the flux subsection
      get_fluxes = True
    elif get_fluxes == True:
      # We are now in the flux subsection.  There may be any number of fluxes
      data = line.split()
      try:
        wavelength = float(data[0][:-2])
        index = 'mm'+str(int(wavelength*10))
        try:
          flux = float(data[6])
          this_source[index] = flux
        except IndexError, details:
          module_logger.debug("Could not parse flux due to index error")
          module_logger.debug(details)
        except ValueError, details:
          # module_logger.error("Warning: no flux in:\n",line,"\nfor",current_source)
          module_logger.debug("Warning: no flux in:\n{}\nfor {}".format(line, current_source))
          module_logger.debug(details)
      except IndexError, details:
        module_logger.debug("Warning: no wavelength in: {",line,"} for",current_source)
        break
      except ValueError, details:
        module_logger.debug("Could not parse wavelength data:",data[0][:-2])
        module_logger.debug("in:",line)
        module_logger.debug(details)
        break
  module_logger.info('{} sources processed'.format(len(cal_data)))
  return cal_data

def make_VLA_pickle_file(cals):
  """
  Put the VLA calibrator dictionary in a pickle file

  @param cals : calibrator data dictionary
  @type  cals : dictionary of dictionaries

  @return: True
  """
  dbfile = open(cal_dir+'VLA_cals','w')
  pickle.dump(cals,dbfile)
  dbfile.close()
  return True

def get_cal_dict():
  """
  Returns the VLA calibrator dictionary
  """
  try:
    dbfile = open(cal_dir+'VLA_cals','r')
  except IOError:
    module_logger.error("file VLA_cals could not be found or else not read")
    # get fresh data
    data = get_VLA_calibrators()
  else:
    data = pickle.load(dbfile)
    dbfile.close()
  return data

def get_cal_data(source):
  """
  Return calibrator data for a source.

  This will handle names like BHHMM+DDd and JHHMM+DDd and 3CNNN.
  If it is an IAU name without a prefix it will try both B and J.
  For now, it doesn't handle other names, like "Virgo A"

  Note
  ====
  This will use the pickle file if there is one.  Otherwise it will get
  the data from NRAO.

  @param source : source name
  @type  source : str

  @return: dictionary
  """
  data = get_cal_dict()
  if source[0].lower() == 'j':
    # This is a IAU J designation
    name = source[1:]
    try:
      source_data = data[name]
    except Exception, details:
      module_logger.debug("Could not find data for J"+name)
      source_data = None
  elif source[0].lower() == 'b':
    # This is an IAU B designator
    try:
      namefile = open(cal_dir+'B_names','r')
      bnames = pickle.load(namefile)
    except IOError, details:
      module_logger.debug("File B_names not found.")
      bnames,cat_3C = VLA_name_xref(data)
    name = source[1:]
    source_data = data[bnames[name]]
  elif source[0:2].lower() == '3c':
    # This is a 3C catalog name
    try:
      namefile = open(cal_dir+'3C_names','r')
      bnames = pickle.load(namefile)
    except IOError, details:
      module_logger.debug("File 3C_names not found.")
      bnames,cat_3C = VLA_name_xref(data)
    source_data = data[cat_3C[source]]
  else:
    # ambiguous
    if data.has_key(source):
      # Matches a J designator
      source_data = data[source]
    elif bnames.has_key(source):
      # Try the B names
      source_data = data[bnames[source]]
      source_data['jname'] = source
    else:
      module_logger.error("No match for",source)
      source_data = None
  return source_data

def VLA_name_xref(cal_data):
  """
  Make cross-reference lookup for IAU B names and 3C catalog numbers.

  Two dictionaries are returned.  The first has IAU B epoch names as keys
  and the second has 3C names as keys.  Both return the corresponding J name.

  @param cal_data : VLA calibrator data dictionary
  @type  cal_data : dictionary of dictionaries

  @return: tuple of dictionaries
  """
  Bname_dict = {}
  cat_3C_dict = {}
  keys = cal_data.keys()
  keys.sort
  for key in keys:
    if cal_data[key].has_key('bname'):
      Bname_dict[cal_data[key]['bname']] = key
    if cal_data[key].has_key('cat3c'):
      cat_3C_dict[cal_data[key]['cat3c']] = key
  return Bname_dict,cat_3C_dict

def make_Bname_pickle_file(Bnames):
  """
  Create a picke file with a Bname to Jname cross-reference

  @param Bnames : cross-reference data

  @return: True
  """
  dbfile = open(cal_dir+'B_names','w')
  pickle.dump(Bnames,dbfile)
  dbfile.close()
  return True

def make_3C_pickle_file(names):
  """
  Create a picke file with a 3C to Jname cross-reference

  @param names : cross-reference data

  @return: True
  """
  dbfile = open(cal_dir+'3C_names','w')
  pickle.dump(names,dbfile)
  dbfile.close()

def IAU_name_parts(name):
  """
  Returns the declination part of an IAU source name
  """
  index = max(name.find('-'),name.find('+'))+1
  ra_part = name[:index-1]
  dec_part = name[index:]
  return ra_part,dec_part

def match_IAU_name(name,name_list):
  """
  Try to match a short IAU name to one in a catalogue

  This is done by rounding the catalogue name as much as necessary.  For
  example, if the name is 1012+53, it is one character short of the coorrect
  name 1012+531. If the name is 1012+5307, it is one character over.
  """
  missing = 8-len(name)
  module_logger.debug("match_IAU_name: missing =",missing)
  if missing == 0:
    # Just see if the name is in the list
    try:
      index = name_list.index(name)
    except ValueError, details:
      module_logger.error("match_IAU_name: could not find",name,"in name list")
      return None
    else:
      return name_list[index]
  # The name is too long or too short
  factor = pow(10.,missing)
  module_logger.debug("match_IAU_name: multiplier =",factor)
  name_ra,name_dec = IAU_name_parts(name)
  module_logger.debug("match_IAU_name: name ra,dec parts>",name_ra,name_dec)
  for candidate in name_list:
    key_ra,key_dec = IAU_name_parts(candidate)
    if key_ra == name_ra:
      module_logger.debug("match_IAU_name: matched ra part to>",candidate)
      if missing > 0:
        # name too small; round the key and see if it is close enough
        if abs(  int(name_dec)
               - int(round(float(key_dec)/factor)) ) < 2:
          return candidate
      elif missing < 0:
        # name too big; round the name and see if it close enough
        if abs(  int(round(float(name_dec)*factor))
               - int(key_dec)                      ) < 2:
          return candidate
  # Went through the whole list without a match
  return None

def fix_IAU_name(name):
  """
  Handles cases where the IAU designator is too long or too short.

  Returns a name suitable for an index into the dictionaries used in
  this module.

  The calling routine needs to ensure that name at least looks like
  an IAU designation.
  """
  # Get the cal data dictionary
  cal_data = get_cal_dict()
  # Find out how to handle the name
  if name[0].lower() == "j" or name[0].lower() == "b":
    # Note which type and truncate the first character
    IAU_type = name[:1].lower()
    name = name[1:]
  else:
    # Leave name unchanged but IAU type is ambiguous
    IAU_type = None
  if len(name) == 8:
    # Right length.  It would be a fluke to have a wrong format
    return name,IAU_type
  else:
    # too short; find the best match.
    # The dec_part could be rounded or truncated
    if IAU_type == "j":
      # Get data keyed to IAU J names
      J_names = get_cal_dict().keys()
      return match_IAU_name(name,J_names),"j"
    elif IAU_type == "b":
      # Get data keyed to IAU B names
      B_names,cat_3C = VLA_name_xref(cal_data)
      return match_IAU_name(name,B_names),"b"
    else:
      # IAU type unknown, try both
      J_names = get_cal_dict().keys()
      result = match_IAU_name(name,J_names)
      if result:
        return result,"j"
      B_names,cat_3C = VLA_name_xref(cal_data)
      result = match_IAU_name(name,B_names)
      if result:
        return result,"b"
      # No match
      return None,None

def Jnames_to_B(Bnames_dict):
  """
  Returns a dictionary for looking up Bnames given Jnames

  @param Bnames_dict : Jnames keyed to Bnames
  @type  Bnames_dict : dictionary

  @return: Bnames keyed to Jnames
  """
  Jnames_dict = {}
  for k, v in Bnames_dict.iteritems():
    Jnames_dict[v] = k
  return Jnames_dict

def get_flux_data(src_data):
  """
  Gets flux as a function of frequency for the source

  @param src_data : data returned from a call to get_cal_data

  @return: (list of freqs, list of fluxes)
  """
  S = {}
  for key in src_data.keys():
    if key[0:2] == "mm":
      if src_data[key] != None and src_data[key] != 0.0:
        mm = int(key[2:])
        f = 300./mm
        S[f] = src_data[key]
  f = S.keys()
  f.sort()
  fluxes = []
  for k in f:
    fluxes.append(S[k])
  return f,fluxes

def interpolate_flux(freqs,fluxes,freq):
  if len(freqs) < 1:
    return None
  elif len(freqs) == 1:
    return fluxes[0]
  elif len(freqs) > 1:
    (ar,br)=polyfit(freqs,fluxes,1)
    return polyval([ar,br],freq)
