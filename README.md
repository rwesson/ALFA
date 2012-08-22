ALFA
====

Automatic line fitting algorithm.

The program reads in an observed spectrum, and a list of lines that might be expected to be present in that spectrum (at the moment, the observed spectrum is a short exposure on NGC 6543, the Cat's Eye Nebula, and the line list is the observed wavelengths of lines measured in a long exposure of the Cat's Eye).

The program attempts to match the observations as follows

1. it generates a set of spectra using the wavelengths of the line list and random fluxes
2. it calculates the RMS of the residuals for each of the synthetic spectra
3. from the spectra with the lowest RMS, it creates new synthetic spectra by averaging random pairs
4. a small fraction of the lines in the new spectra have their flux altered
5. steps 2-4 are repeated many times, ideally resulting in increasingly good fits to the spectrum being generated
