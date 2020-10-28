# LGCATuringPatterns
Repo for Julia scripts needed for producing the results of the LGCA Turing pattern project


StateSpaceProducer.jl - Produces the asymptotic state space for the 21 topologies. 

SimulationFunctions.jl - Contains the functionality to run the simulations of a given map. 

PowerSpectrum.jl - Contains the functionality to produce the power spectrum data for a given map. Dependency on CoreSimulator.jl 

SteadyStateFinder.jl - Uses numerical methods to find the steady states of the mean field approximated system of a given map.  

LinearStabilityAnalysis.jl - Provides the dispersion relation of a given map using the steady states found using SteadyStateFinder.jl

GridPlotter.jl - A script that plots the 3 trial simulations, the dispersion relation (using the first steady state if more than one) and the powerspectrum (using 100 trials) for a given map. 
