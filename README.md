ALFA
====

Automatic line fitting algorithm.

ALFA can identify and fit hundreds of lines in emission line spectra in just a few seconds.  It does this using a genetic algorithm to optimise the line parameters, which works as follows:

1. A population of synthetic spectra is generated using a reference line catalogue.
2. The goodness of fit for each synthetic spectrum is calculated.  The best sets of parameters are retained and the rest discarded.
3. A new population of synthetic spectra is obtained by averaging pairs of the best performers.
4. A small fraction of the parameters of the lines in the new generation are randomly altered.
5. The process repeats until a good fit is obtained.

The code is currently optimised for optical spectra but the algorithm is entirely indifferent to the wavelength range and resolution of the spectra to be analysed.  The only requirement to get a good fit out is to create a meaningful reference line catalogue.

Installation
============

To install ALFA, unpack the source files somewhere convenient, then type

  make
  sudo make install

If you don't have the necessary permissions to execute the second command, then change the directory definitions in the make file, edit source/alfa.f95 and source/alfa_cube.f95 to change the values of skylinelistfile, stronglinelistfile and deeplinelistfile, and then type

  make install
