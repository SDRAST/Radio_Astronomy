Radio\_Astronomy
================

Radio astronomical engineering support, radio source and other calibrators.

The code can be cloned from `this site <https://github.com/SDRAST/Radio_Astronomy>`_

The `documentation <https://sdrast.github.io/Radio_Astronomy/>`_
was generated with Sphinx. 

Examples
--------
**Antenna Temperature (K) of Venus Today for a 34-m Antenna at 8.4 GHz**

The radio flux of Venus depends on its distance and the frequency. The antenna
temperature registered by the radio telescope depends on its size and efficiency::

    In [1]: from math import pi
    In [2]: from datetime import datetime
    In [3]: from Radio_Astronomy import antenna_gain
    In [4]: from Radio_Astronomy.radio_flux import get_planet_flux
    In [5]: antenna_gain(0.7, pi*(34/2)**2) * get_planet_flux('Venus', 8.4, datetime.now())
    Out[5]: 19.696163420934074

