
# System identification exapmle

This is an example of a simple mechanical system identification. The system
consist on a Mass-Spring-Damper system.

A simulation is done for a time span with a fixed sample period.
Position, velocity and acceleration is recorded.

The simulation is repeated with a periodic force input to introduce some
fixed frequency content.

Finally, the system is excited with an input signal generated with an
identification toolbox. It contains all the frequencies in a span.
The FRF of the system is calculated and displayed.

# Toolbox used

The toolbox created at [MECO][meco-group] group of KU Leuven.

# Author

Sergio Portoles

# Debug

Some problems were encountered during the FFT of the input signal.
Maybe the rest of the signals have the problem too.

[meco-group]: https://www.mech.kuleuven.be/en/pma/research/meco "MECO at KU Leuven"
