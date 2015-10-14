# -*- coding: utf-8 -*-
"""
Data for 86.1 GHz from
Ulich, B. L.; Davis, J. H.; Rhodes, P. J.; Hollis, J. M.
IEEE Trans. Ant. and Prop., AP-28, 367-377 (1980).

Data for 14.5 GHz from
Baars, J.W.M., Mezger, P.G., Wendker, H.
Zeitschrift fuer Astrophysik, 61, 134-143 (1965)
"""

from math import exp, pow, log10
import ephem
import numpy as NP
from scipy.optimize import leastsq

from Radio_Astronomy import flux

# planets recognized by module ephem
Planets = ['Jupiter', 'Mars', 'Mercury', 'Moon', 'Neptune', 'Pluto',
           'Saturn', 'Sun', 'Uranus', 'Venus']
diag = False

params = {"Virgo": [6.541, -1.289],
          "Omega": [4.056, -0.378],
          "Orion": [3.317, -0.204]}

def radio_flux(source,freq):
  """
  Flux density in Jy of standard calibrators.

  From Radio Star Flux Density Expressions for Accurate Antenna Gain
  Measurements by E.P. Ekelman, COMSAT Laboratories,
  0-7803-5639-W99/  IEEE 1999.

  @param source : name of the source ('Virgo', 'Omega', or 'Orion')
  @type  source : str

  @param freq : frequency in GHz
  @param freq : float

  @return: the flux density in Jy.
  """
  p1 = params[source][0]
  p2 = params[source][1]
  return pow(10,(p1+p2*log10(1000*freq)))

def subearth_longitude(pl):
  """
  Longitude of planet facing earth from heliocentric longitude

  This uses the elongation to compute the angle between the Sun and the
  Earth as seen from the planet and then subtracts or adds that to the
  heliocentric longitude to get the sub-earth longitude.

  @param pl : planet
  @type  pl : ephem.Planet() instance
  """
  return pl.hlong+pl.elong*(pi/2/abs(elong) - 1)

def log_cubic(x,p):
  """
  2-order polynomial in log(x)
  """
  logx = log10(x)
  return p[0] + p[1]*logx + p[2]*logx**2 + p[3]*logx**3

def err_func(p,x,y,yerr,function):
  """
  Difference between data and a model
  """
  #print "p is",type(p),"of length",len(p)
  #print "x is",type(x),"of length",len(x)
  #print "y is",type(y),"of length",len(y)
  fit = function(x,p)
  #print "fit is",type(fit),"of length",len(fit)
  return (y - function(x,p))/yerr**2
  
def planet_brightness(planet, freq):
  """
  Brightness temperature of a planet

  @param planet : name of the planet (or 'Sun')
  @type  planet : str

  @param freq : frequency (or frequencies) in GHz
  @type  freq : float or array of float

  @return: float or array of float
  """
  if planet == 'Venus':
    # Butler et al, Icarus, 154, 226B (2001) for 4.86-22.46 GHz
    # Ulich et al., IEEE Proc. Ant. Prop., AP-28, 367-377 (1980) for 86.1 GHz
    # Yefanov et al. Radiofizika, 13, 219 (1970) for 37.5 and 138.9
    # Baars et al., Z.f.Astroph. 61, 134 (1965) for 14.5 GHz
    # 
    # This fits the data down to 4 GHz but not below
    freqs  = NP.array([ 22.46, 14.94,  8.44,  4.86, 14.5,  86.1,  37.5, 138.9])
    Tb     = NP.array([505.2, 565.9, 657.5, 679.9, 480.0, 357.5, 495.,  290. ])
    #sig_Tb= array([ 25.3,  17.0,  13.2,  13.6,  50.0,  13.1,  20.0,  25.0])
    sig_Tb= NP.array([ 25.3,  17.0,  13.2,  13.6,  100.0,  13.1, 100.0, 100.0])
    pinit = [750., -100., 0., 0.]
    out = leastsq(err_func, pinit,
                  args=(freqs,Tb,sig_Tb,log_cubic), full_output=1)
    Tb_fit = log_cubic(freq,out[0])
    mask = freq <= 4.0
    if type(mask) == bool:
      if freq <= 4.0:
        return 686.0
      elif freq > 75.0:
        return 351.0
      else:
        return Tb_fit
    else:
      if mask.any():
        Tb_fit = Tb_fit*(freq > 4.0) + mask*686.0
      mask = freq >= 75.0
      if mask.any():
        Tb_fit = Tb_fit*(freq < 75.0) + mask*351.0
      return Tb_fit
  elif planet == 'Sun':
    # Ulich et al., IEEE Proc. Ant. Prop., AP-28, 367-377 (1980) for 86.1 GHz
    freqs  = [ 86.1]
    Tb     = [7914.]
    sig_Tb = [ 192.]
    pass
  elif planet == 'Jupiter':
    freqs  = [ 14.5,  86.1]
    Tb     = [157.0, 179.4]
    sig_Tb = [ 14.0,   4.7]
    pass
  elif planet == 'Saturn':
    freqs  = [86.1]
    Tb     = [153.4]
    sig_Tb = [4.8]
    pass
  
def get_planet_flux(planet,freq,date):
  """
  Flux of planet calibrators

  Sources
  =======
  The Flux Density of the Strongest Thermal Radio Sources at 14.5 GHz
  by J.W.M. Baars, P.G. Mezger and H. Wendker

  @param planet : planet name
  @type  planet : str

  @param freq : frequency in GHz
  @type  freq : float

  @param date : date and time of observation
  @type  date : datetime.datetime() instance

  @return: flux in Jy
  """
  source = planet.capitalize()
  if source == 'Venus':
    pl = ephem.Venus()
    Tb = planet_brightness('Venus', freq)
    sig_Tb = 5
    ref_Tb = 'model'
  elif source == 'Jupiter':
    pl = ephem.Jupiter()
    Tb = 157
    sig_Tb = 14
    ref_Tb = 'Baars1965'
  elif source == 'Saturn':
    pl = ephem.Saturn()
    Tb = 0.94*157
    sig_Tb = 0.94*14
    ref_Tb = 'Welch1966'
  elif source == 'Sun':
    pl = ephem.Sun()
    Tb = 5800
    sig_Tb = 100
    ref_Tb = ''
  elif source == 'Mars':
    # 
    pl = ephem.Mars()
    pl.compute(date)
    Tb = 190
    sig_Tb = 12
    # This corrects for the insolation as a function of distance from the Sun
    R = pl.sun_distance
    Tb *= (1.524/R)**2
    sig_Tb *= (1.524/R)**2
    ref_Tb = 'Dent1965'
  elif source == 'Mercury':
    pl = ephem.Mercury()
    pl.compute(date)
    D = pl.phase/100
    l = 300/freq # mm
    Tb = 330 * pow(10,(0.1-0.4*D)/l)
    ref_Tb = 'Klein1970'
  if pl.name != 'Mercury' and pl.name != 'Mars':
    pl.compute(date)
  diameter = 2*pl.radius           # angle in radians
  if diag:
    print "Tb =",Tb, ", radius =", radius
  # This converts brightness temperature and size to flux.
  return flux(Tb, freq, diameter)

def get_calibrator_flux(source, freq, date):
  """
  Get the flux of a calibrator by name, frequency and date.

  This first tries to see if a planet matches and if so computes the flux
  for that.  It then tries a 3C catalogue name.  If there is a match, it
  tries to find it in the Michigan catalogue and interpolates or extrapolates
  the flux from that.  Otherwise, it tries to find it in the VLA calibrator
  catalogue.

  It does not yet handle 'B' or 'J-\' IAU names, though that would be a
  trivial extension.  The Michigan catalogue uses 'B' names.

  @param source : name of calibrator
  @type  source : str

  @param freq : frequency in GHz
  @type  freq : float

  @param date : date for flux estimate, for variable sources
  @type  date : datetime.datetime() instance

  @return: flyx in Jy (float) and origin of flux data (str)
  """
  if diag:
    print "Processing",source,"for",freq,"GHz at",date.ctime()
  try:
    # planet?
    Planets.index(source)
    flux = get_planet_flux(source, freq, date)
    ref = "Planet"
  except ValueError:
    # not a planet
    if source[1] == 'C' or source[0] == 'J' or source[0] == 'B':
      calibrator = Ephem.Quasar(source)
      if diag:
        print source,"=",calibrator.Jname,"=",calibrator.Bname
      flux = calibrator.get_flux(freq,date)
      if diag:
        print ref,"flux is",flux
    else:
      try:
        # Maybe its a J name without the J
        calibrator = Ephem.Quasar(source)
        flux = calibrator.get_flux(freq,date)
        if diag:
          print ref,"flux is",flux
      except:
        # Need to handle more cases
        if diag:
          print "Could not handle",source
        flux = None
        ref = None
  return flux, ref

def galactic_BG(f):
  """
  Average galactic background temperature

  Zarka et al JGR, 109, A09S15 (2004)

  @param f : frequency in MHz
  @type  f : float
  """
  def tau(f):
    """
    Optical depth of the Galaxy

    @param f : frequency in MHz
    @type  f : float
    """
    return 5.0*f**(-2.1)

  Ig = 2.48e-20
  Ieg = 1.06e-20
  return Ig*f**(-0.52)*( (1-exp(-tau(f)))/tau(f) ) + Ieg*f**(-0.8)*exp(-tau(f))

if __name__ == "__main__":
  from pylab import *
  freqs  = NP.array([ 22.46, 14.94,  8.44,  4.86, 14.5,  86.1,  37.5, 138.9])
  Tb     = NP.array([505.2, 565.9, 657.5, 679.9, 480.0, 357.5, 495.,  290. ])
  semilogx(freqs,Tb,'.')
  x = linspace(1,100,100)
  y = planet_brightness("Venus", x)
  semilogx(x,y,'-')
  grid()
  show()
