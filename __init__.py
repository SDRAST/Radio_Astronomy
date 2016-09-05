# -*- coding: utf-8 -*-
"""
Radio_astronomy - basic radio astronomy functions

Function Groups
===============

Antenna
-------
Functions for antenna properties::
   antenna_gain(aperture_efficiency,geometrical_area)
   antenna_solid_angle(aperture_efficiency,geometrical_area,wavelength)
   antenna_temperature(flux,effective_area)
   beam_efficiency(antenna_solid_angle, beam_solid_angle)
   beam_solid_angle(*args)
   directivity(aperture_efficiency, geometrical_area, wavelength)
   HPBW(feedtaper,wavelength,diameter)
   ruze_loss_factor(surface_rms,wavelength)

Noise
-----
Functions related noise temperature::
   antenna_temp(antenna_gain,flux_density)
   noise_figure(Tsys)
   noise_power(Tsys,bandwidth)
   rms_noise(Tsys,bandwidth,integration_time)
   SNR(Ta,Tsys,bandwidth,integration_time)

Radio Sources
-------------
Relating flux, antenna temperature, etc.::
   flux(Tb,freq,angular_diameter)
   janskyPQ(Jy)

Conversions
-----------
Converting between units, log to linear, etc.::
   dB(gain)
   dBm(power)
   dbm_to_dbuv(dbm)
   dbuv_to_uV(dbuv)
   dbm_to_v(dbm)
   dBm_to_watts(dbm)
   dbuv_to_dbm(dbuv)
   dbuv_to_dmw_per_sq_m(dbuv)
   dbuv_to_v(dbuv)
   dmw_per_sq_m_to_dbuv(dbmw)
   gain(dB)
   v_to_dbuv(v)
   volts_to_watts(V)
   v_to_dbm(V)
   watts_to_volts(W)

Miscellaneous
-------------
Left-overs::
   freq_to_chan(frequency,bandwidth,n_chans)

"""
import math
try:
  import Physics
  from Physics import c, wavenumber
  from Physics.Radiation.Continuum import BB_intensity
  
  def janskyPQ(Jy):
    """
    jansky unit (Jy) as a PhysicalQuantity() class instance

    @param Jy : quantity expressed in Jy
    @type  Jy : float

    @return: PhysicalQuantity instance
    """
    return Physics.pq(Jy*1e-26,"W") \
          /Physics.pq(1,'m')/Physics.pq(1,'m') \
          /Physics.pq(1,'Hz')

  Jy = janskyPQ(1).getValue()
except ImportError:
  Jy = 1e-26
  c = 299792458.0
  class Physics:
    c = 299792458.0
    h = 6.6260755e-34
    k = 1.380658e-23
  def wavenumber(wvln):
    return 1./wvln
  def BB_intensity(T,f):
    left_term = 2*h*math.pow(f,3)/math.pow(c,2)
    exponent = h*f/k/T
    if exponent < 700:
      right_term = math.pow(math.e,exponent) - 1
      result = left_term/right_term
    else:
      result = 1e-304
    return result

import numpy
import sys
from scipy.special import sici

euler = 0.5772156649 # Euler constant
python_version = \
  "python"+str(sys.version_info.major)+"."+str(sys.version_info.minor)

cal_dir = "/usr/local/lib/"+python_version+"/DSN-Sci-packages/Radio_Astronomy/"

# -------------------------- classes -------------------------------

class Dipole():
  """
  Dipole antenna
  """
  def __init__(self,length,radius):
    """
    Create a dipole antenna instance

    @param length : dipole length (m)
    @type  length : (numpy array of) float
    
    @param radius : radius of dipole element (m)
    @type  radius : (numpy array of) float
    """
    self.length = length
    self.radius = radius

  def impedance(self, wvln, Z):
    """
    Impedance of a dipole antenna
    
    Reference::
      https://en.wikipedia.org/wiki/Dipole_antenna#General_impedance_formulas
    
    @param wvln : wavelength (m)
    @type  wvln : (numpy array of) float
    
    @param Z : impedance of the medium, 376.73 ohm for free space
    """
    k = 2*math.pi*wavenumber(wvln)
    x = k*self.length
    Si,Ci = sici(x)
    SiTwo,CiTwo = sici(2*x)
    Sia, Cia = sici(2*k*self.radius*self.radius/self.length)
    multiplier = Z/(2*math.pi*numpy.sin(x/2.)**2)

    R = multiplier * (euler + numpy.log(x) - Ci
                      +0.5*numpy.sin(x)*(SiTwo - 2*Si)
                      +0.5*numpy.cos(x)*(euler + numpy.log(x/2.) + CiTwo -2*Ci))
    X = (multiplier/2.) * (2*Si + numpy.cos(x)*(2*Si - SiTwo)
                         - numpy.sin(x)*(2*Ci - CiTwo - Cia))
    return R, X



# ------------------------ global methods --------------------------

def antenna_gain(aperture_efficiency,geometrical_area):
  """
  Antenna gain in K/Jy.

  Given the geometrical aperture area in m^2 and the aperture
  efficiency, this returns the gain in K/Jy

  @param aperture_efficiency :
  @type  aperture_efficiency : float

  @param geometrical_area : typically pi*(diam/2)^2, same units as wavelength
  @type  geometrical_area : float

  @return: float
  """
  return aperture_efficiency*geometrical_area/2/Physics.k/1e26

def antenna_solid_angle(aperture_efficiency,geometrical_area,wavelength):
  """
  Antenna solid angle of an antenna, assuming radiation efficiency = 1

  Given the aperture efficiency (dimensionless) and the geometrical
  area and wavelength in the same units (say, meters), this returns the
  beam solid angle in steradians.

  If the antenna has some losses, like not being totally reflective,
  then the radiation efficiency is not 1 and it should be multiplied
  with this result.

  @param aperture_efficiency :
  @type  aperture_efficiency : float

  @param geometrical_area : typically pi*(diam/2)^2, same units as wavelength
  @type  geometrical_area : float

  @param wavelength : same units as antenna  diameter
  @type  wavelength : float

  @return: float
  """
  return math.pow(wavelength,2)/geometrical_area/aperture_efficiency

def antenna_temp(antenna_gain,flux_density):
  """
  Antenna temperature of a source from flux and antenna gain.
  
  Given the gain of the antenna (K/J), and the flux density of a
  source in Jy, this returns the antenna temperature of the source.

  @param antenna_gain : antenna gain in K/Jy
  @type  antenna_gain : float
  
  @param flux_density : flux density in Jy
  @type  flux_density : float

  @return: antenna temperature in K (float)
  """
  return antenna_gain*flux_density

def antenna_temperature(flux,effective_area):
  """
  Convert flux to antenna temperature

  @param flux : source flux in Jy
  @type  flux : float

  @param effective_area : in m^2
  @type  effective_area : float

  @return: antenna temperature in K
  """
  return flux*Jy*effective_area/Physics.k/2

def beam_efficiency(antenna_solid_angle, beam_solid_angle):
  """
  Fraction of flux in the main beam that registers at the receiver.

  This is generally larger than the aperture efficiency by a factor
  of 1.2 - 1.4.
  """
  return beam_solid_angle/antenna_solid_angle

def beam_solid_angle(*args):
  """
  Solid angle of the main beam out to the first null.

  This is an approximation assuming a Gaussian beamshape.

  This is also an alias for antenna_solid_angle for backwards compatibility.
  The number of arguments determines which it is.
  """
  if len(args) == 2:
    HPBWx, HPBWy = args
    return 1.13*HPBWx*HPBWy
  elif len(args) == 3:
    aperture_efficiency,geometrical_area,wavelength = args
    return antenna_solid_angle(aperture_efficiency,geometrical_area,wavelength)
  else:
    return None

def dB(gain):
  """
  Converts a linear gain into dB

  @param gain : power out / power in ratio
  @type  gain : float

  @return: dB (float)
  """
  if type(gain) == float:
    return 10*math.log10(gain)
  elif type(gain) == list:
    return 10*numpy.log10(numpy.array(gain))
  else:
    return 10*numpy.log10(gain)

def dBm(power):
  """
  Converts a power in W to db(milliwatts)

  @param power : power in W
  @type  power : float

  @return: dBm (float)
  """
  if type(power) == float or type(power) == int:
    return 10*math.log10(power*1000.)
  elif type(power) == numpy.ndarray:
    return 10*numpy.log10(power*1000.)
  else:
    return 10*numpy.log10(numpy.array(power)*1000.)
    
def dbm_to_dbuv(dbm):
  """
  Convert dB(milliwatts) to dB(microvolts)

  This is the same as dbm + 106.98 except pedagogically clearer

  @param dbm : power in dBm
  @type  dbm : float

  @return: dbuv (float)
  """
  W = dBm_to_watts(dbm)
  uV = watts_to_volts(W)*1.e6
  if type(uV) == float or type(uV) == int:
    return 20.*math.log10(uV)
  else:
    return 20.*numpy.log10(uV)

def dbm_to_v(dbm):
  """
  Convert dBm to V across a 20-ohm load

  @param dbm : power in dBm
  @type  dbm : float

  @return: V (float)
  """
  return numpy.sqrt(10**(dbm/10.)/1000*50)

def dBm_to_watts(dbm):
  """
  Convert dBm into watts

  @param dbm : power in dBm
  @type  dbm : float

  @return: watts (float)
  """
  if type(dbm) == float or type(dbm) == int:
    return math.pow(10.,dbm/10.)/1000.
  else:
    return numpy.power(10.,dbm/10.)/1000.


def dbuv_to_dbm(dbuv):
  """
  Convert dBuv to dBm.
  
  This is the same as dbm - 106.98 except pedagogically clearer

  @param dbuv : dB(uV) as float
  """
  V = math.pow(10.,dbuv/20.)/1.e6
  W = volts_to_watts(V)
  return dBm(W)

def dbuv_to_dmw_per_sq_m(dbuv):
  """
  Convert dB(uV) to dB(mW/m^2)
  
  from http://www.ahsystems.com/notes/RFconversions.php:
  dBmW/m2 = dBmV/m - 115.8

  @param dbuv : dB(uV)
  @type  dbuv : float

  @return: dBm/m^2 as float
  """
  return dbuv - 115.8

def dbuv_to_uV(dbuv):
  """
  dBuv (dB microvolts) to V (volts)

  @param dbuv : dB(uV)
  @type  dbuv : float

  @return: uV (float)
  """
  return (10**(dbuv/20.))

def dbuv_to_v(dbuv):
  """
  dBuv (dB microvolts) to V (volts)

  @param dbuv : dB(uV)
  @type  dbuv : float

  @return: V (float)
  """
  return dbuv_to_uV(dbuv)/1e6

def dmw_per_sq_m_to_dbuv(dbmw):
  """
  Convert dB(mW/m^2) to dB(uV)
  
  from http://www.ahsystems.com/notes/RFconversions.php:
  dBmW/m2 = dBmV/m - 115.8

  @param dbmw : dB(mW/m^2)
  @type  dbmw : float

  @return: dB(uV) as float
  """
  return dbmw + 115.8

def directivity(aperture_efficiency, geometrical_area, wavelength):
  """
  Forward gain of an antenna over an isotropic antenna as a factor

  Given the aperture efficiency (dimensionless) and the geometrical
  area and wavelength in the same units (say, meters), this returns
  the forward or directive gain.

  @param aperture_efficiency :
  @type  aperture_efficiency : float

  @param geometrical_area : typically pi*(diam/2)^2, same units as wavelength
  @type  geometrical_area : float

  @param wavelength : same units as antenna  diameter
  @type  wavelength : float

  @return: float
  """
  return 4*math.pi * aperture_efficiency \
                   * geometrical_area \
                   / math.pow(wavelength,2)

def forward_gain(aperture_efficiency, geometrical_area, wavelength):
  """
  Directivity of an antenna in dB.

  See 'directivity' for parameters.
  """
  return dB(directivity(aperture_efficiency, geometrical_area, wavelength))
  
def flux(Tb,freq,angular_diameter):
  """
  Flux received from a source

  @param Tb : brightness temperature in K
  @type  Tb : float

  @param freq : frequency in GHz
  @type  freq : float

  @param angular_diameter : source diameter in radians
  @type  angular_diameter : float

  @return: flux in Jy
  """
  I = BB_intensity(Tb,freq*1e9)
  solid_angle = math.pi*(angular_diameter/2)**2 # small angle approximation
  return I*solid_angle/Jy

def freq_to_chan(frequency,bandwidth,n_chans):
    """
    Returns the channel number where a given frequency is to be found.

    @param frequency : same units as bandwidth
    @type  frequency : float
    
    @param bandwidth : same units as frequency
    @type  bandwidth : float
    
    @param n_chans : number of channels in the band
    @type  n_chans : int
    """
    if frequency < 0:
        frequency = bandwidth + frequency
    if frequency > bandwidth:
      raise RuntimeError("that frequency is too high.")
    return round(float(frequency)/bandwidth*n_chans) % n_chans

def gain(dB):
  """
  Convert dB into gain

  @param dB : decibels
  @type  dB : float

  @return: ratio (float)
  """
  return 10.**(dB/10.)

def HPBW(feedtaper,wavelength,diameter):
  """
  Half-power beamwidth estimate

  @param feedtaper : feed pattern amplitude at edge of primary, in dB
  @type  feedtaper : float

  @param wavelength : in same units as diameter
  @type  wavelength : float

  @param diameter : of main aperture, in same units as wavelength
  @type  diameter : float

  @return: HPBW in radians (float)
  """
  return (1.02 + 0.0135 * feedtaper) * wavelength/diameter

def noise_figure(Tsys):
  """
  Returns the noise figure in dB given a system temperature in K

  @param Tsys : system temperature (K)
  @type  Tsys : float

  @return: noise figure (excess over ambient) in dB
  """
  return 10*math.log10(1 + Tsys/290.)

def noise_power(Tsys,bandwidth):
  """
  Returns the noise power in W

  @param Tsys : system temperature in K
  @type  Tsys : float
  
  @param bandwidth : bandwidth in Hz
  @type  bandwidth : float

  @return: noise power in W (float)
  """
  return Physics.k*Tsys*bandwidth

def rms_noise(Tsys,bandwidth,integration_time):
  """
  Receiver r.m.s. noise temperature.

  Given the system temperature in K, the bandwidth in Hz, and the
  integration time in sec, returns the r.m.s. noise of the sample.

  @param Tsys : system temperature in K
  @type  Tsys : float

  @param bandwidth : bandwidth in Hz
  @type  bandwidth : float

  @param integration_time : integration time in sec
  @type  integration_time : float

  @return: float
  """
  return Tsys/math.sqrt(bandwidth*integration_time)

def ruze_loss_factor(surface_rms,wavelength):
  """
  Loss due to reflecting surface irregulaties in a paraboloid
  
  Given a surface r.m.s. deviation from perfect, this returns the
  efficiency of a reflector relative to one with a perfect surface at
  the specified wavelength.  surface_rms and wavelength must be in the
  same units.

  The aperture (or beam) efficiency is this times the ideal aperture
  (or beam) efficiency considering blockage and beam taper.

  @param surface_rms : r.m.s. of the surface deviation from a parabola
  @type  surface_rms : float

  @param wavelength : in the same units as surface_rms.
  @type  wavelength : float

  @return: 0 < float < 1
  """
  return math.exp(-(4*math.pi*surface_rms/wavelength)**2)

def SNR(Ta, Tsys, bandwidth, integration_time):
  """
  Signal-to-noise ratio for a signal Ta relative to the system noise

  @param Ta : signal antenna temperature (K)
  @type  Ta : float

  @param Tsys : system temperature (K)
  @type  Tsys : float

  @param bandwidth : system bandwidth in Hz
  @type  bandwidth : float

  @param integration_time : integration time
  @type  integration_time : float

  @return: SNR (float)
  """
  return Ta/rms_noise(Tsys,bandwidth,integration_time)

def v_to_dbuv(v):
  """
  Volts to dBuv (dB microvolts)

  @param v : volts
  @type  v : float

  @return: dB(uV)
  """
  return 20*numpy.log10(v*1e6)

def volts_to_watts(V):
  """
  Volts to watts for a 50 ohm load

  @param V : volts
  @type  V : float

  @return: watts
  """
  return V**2/50.

def v_to_dbm(V):
  """
  Volts to dB(mW) for a 50 ohm load

  @param V : volts
  @type  V : float

  @return: dBm
  """
  return dBm(volts_to_watts(V))

def watts_to_volts(W):
  """
  Watts to volts for a 50 ohm load

  @param W : power in watts
  @type  W : float

  @return: watts (float)
  """
  if type(W) == float or type(W) == int:
    return math.sqrt(W*50.)
  else:
    return numpy.sqrt(W*50.)

def delta_f(delta_v,frequency):
  """
  Frequency shift corresponding to a velocity shift.

  For example, 1 km/s = 74 kHz at 22 GHz water line and 6 kHz at the OH line::
   In [1]: from Radio_Astronomy import *
   
   In [2]: delta_f(1e3,22.235e9)/1e3
   Out[2]: 74.167976567309097
   
   In [3]: delta_f(1e3,1776e6)/1e3
   Out[3]: 5.9240983307191799

  @param delta_v : velocity in m/s

  @param frequency : rest frequency in Hz

  @return: frequency shift in Hz
  """
  return (delta_v/Physics.c)*frequency
=======
"""
Pseudo front end for IFs coming from wherever
"""
import logging

from MonitorControl import Port, Beam, ComplexSignal, IF
from MonitorControl.FrontEnds import FrontEnd

module_logger = logging.getLogger(__name__)

class SignalSources(FrontEnd):
  """
  Class for describing the IFs entering the DTO SamGen

  For now this is specific to the signal sources in the PSDG lab.  It will need
  to be generalized.
  """
  def __init__(self, name, inputs=None, output_names=None):
    """
    Initiate a group of input signals

    Generally this is created without inuts.

    @param name : unique name for the signal group
    @type  name : str

    @param inputs : input ports of the signal source group
    @type  inputs : dict Port instance

    @param output_names : name for the output port in a list
    @type  output_names : list of str
    """
    mylogger = logging.getLogger(module_logger.name+".SignalSources")
    FrontEnd.__init__(self, name, inputs=inputs, output_names=output_names)
    self.logger = mylogger
    self.data['frequency'] = 0.320 # GHz
    self.data['bandwidth'] = 0.640
    self.channel={}
    for key in self.outputs.keys():
      index = output_names.index(key)
      self.channel[key] = self.Channel(self,key,
                                        output_names=[key],
                                        signal=Beam('S'+str(index)))
    self.logger.debug(" %s outputs: %s", self, str(self.outputs))

  class Channel(FrontEnd.Channel):
    """
    Class for the path of one IF
    """
    def __init__(self, parent, name, inputs=None, output_names=None,
                 signal=None, active=True):
      """
      Create a specific FrontEnd subclass instance

      @param name : unique name for the signal group
      @type  name : str

      @param inputs : input ports of the signal source group
      @type  inputs : dict Port instance

      @param output_names : name for the output port in a list
      @type  output_names : list of str

      @param signal : the signal handled by this channel
      @type  signal : a Signal class instance
      """
      self.logger = logging.getLogger(parent.logger.name+".Channel")
      FrontEnd.Channel.__init__(self, parent, name, inputs=inputs,
                                  output_names=output_names, active=active)
      self.logger.debug(" initializing SignalSources channel %s", self)
      self.logger.debug(" %s inputs: %s", self, str(inputs))
      self.outputs[name] = Port(self, name, source=None,
                                signal=IF(ComplexSignal(signal,'H'),'I'))
      self.outputs[name].signal['beam'] = signal.name
      self.outputs[name].signal['frequency'] = parent['frequency']*1e9 # Hz
      self.outputs[name].signal['bandwidth'] = parent['bandwidth']*1e9 # Hz
      parent.outputs[name] = self.outputs[name]

class DSNfeSX(FrontEnd):
  """
  The DSN S/X receiver
  """
  def __init__(self, name, inputs=None, output_names=None):
    """
    """
    self.name = name
    mylogger = logging.getLogger(module_logger.name+".SX_fe")
    mylogger.debug(" initializing %s", self)
    mylogger.debug(" %s input channels: %s", self, str(inputs))
    mylogger.debug(" output names: %s", output_names)
    FrontEnd.__init__(self, name, inputs=inputs, output_names=output_names)
    self.logger = mylogger
    self.channel = {}
    keys = self.inputs.keys()
    keys.sort()
    for feed in keys:
      index = keys.index(feed)
      beam_signal = Beam(feed)
      for prop in self.data.keys():
        beam_signal.data[prop] = self.data[prop]
      #beam_signal.name = feed
      self.channel[feed] = self.Channel(self, feed,
                                        inputs={feed: self.inputs[feed]},
                                        output_names=output_names[index],
                                        signal=beam_signal)
    self.logger.debug("%s output channels: %s", self, str(self.outputs))

