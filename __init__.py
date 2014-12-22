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
