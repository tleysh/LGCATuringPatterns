#Power spectrum module
#Dependency on Simulation Functions
include("SimulationFunctions.jl")
using FFTW

m = [1,7]; L = 101; reactionprob = 0.9; totaltime = 250;
MapFormat = SimulationFunctions.TopologyAndMap(9,300);
WholeArray = SimulationFunctions.ThePossibleChannelArrangements()
Intial = SimulationFunctions.InitialArray(2,101)


#Producing the Simulation data for a given map
function ProduceSimDataForGivenMap(MapFormat::Array{Int8,2},m::Array{Int64,1},L::Int64,WholeArray::Array{Array{Array{Int8,2},1},1},reactionprob::Float64,totaltime::Int64)
	SecondTimePoint = SimulationFunctions.InitialArray(2,101)
	FirstTimePoint = SimulationFunctions.Onestep8bit(Intial,m,MapFormat,L,reactionprob,reactionprob,WholeArray)
	for timestep = 1:(totaltime-1)
		FirstTimePoint = SimulationFunctions.Onestep8bit(SecondTimePoint,m,MapFormat,L,reactionprob,reactionprob,WholeArray)
		SecondTimePoint = SimulationFunctions.Onestep8bit(FirstTimePoint,m,MapFormat,L,reactionprob,reactionprob,WholeArray)
	end
	CombinedFirst = SimulationFunctions.Combining8Bitchannels(FirstTimePoint,L)
	CombinedSecond = SimulationFunctions.Combining8Bitchannels(SecondTimePoint,L)
	return CombinedFirst,CombinedSecond
end

#apply a function and map to do the power spectrum
function Power(OneFreq::Float64)
	power = (1/L)*(abs(OneFreq)^2)
	return power
end

#Producing power spectrum data for a single map
function PowerSpectsForEachIndividualSim(MapFormat::Array{Int8,2},m::Array{Int64,1},L::Int64,WholeArray::Array{Array{Array{Int8,2},1},1},reactionprob::Float64,;totaltime = 501,trials = 100)
	AllPowerSpects = Array{Array{Float64,1},1}(undef,trials)
	for i = 1:trials
			SimDataLastTwoSteps = ProduceSimDataForGivenMap(MapFormat,m,L,WholeArray,reactionprob,totaltime)
			#Element wise summation of the last two time steps
			AverageSimData = sum(SimDataLastTwoSteps)/2
			#Define the average
	     AverageData = sum(AverageSimData)/L
		 #Element wise subtraction of average data
	     NormaliseData = AverageSimData.-AverageData
		 #Taking the fft of the data
	     FftOfData = fft(NormaliseData)
		 #Apple map to get power spectrum
	     powerspec = map(OneFreq -> Power(OneFreq),FftOfData)
	     #Power Spectrum halved
	     HalfOfPower = powerspec[1:Int(round(L/2))]
			if sum(HalfOfPower) != 0
	         NormedHalf = HalfOfPower/sum(HalfOfPower)
	          AllPowerSpects[i] = NormedHalf
			else
				NormedHalf = ones(Int(round(L/2)))*(1/(Int(round(L/2))))
				 AllPowerSpects[i] = NormedHalf
			end

	end
    return AllPowerSpects
end


#Producing normalisated Powerspectrum for n number of trials
function Normalisation(PowerSpect::Array{Array{Float64,1},1},trials::Int64,L::Int64)

   AllFrequencies = zeros(Int(round(L/2)))
   for currentfreq = 1:Int(round(L/2))
      onefreq = 0
      for i = 1:trials
         onefreq = onefreq + powerspectsunnorme[i][currentfreq]
      end
      AllFrequencies[currentfreq] = onefreq
   end
   normalisedfreqs = AllFrequencies./trials
   return normalisedfreqs
end


#Alltogether function to go from map to normalised power spectrum
function NormalisedPowerSpectrum(MapFormat::Array{Int8,2},m::Array{Int8,1},L::Int64,WholeArray::Array{Array{Array{Int8,2},1},1},reactionprob::Int64,;totaltime = 501,trials = 100)

	AllPowerFirstSpecies = PowerSpectsForEachIndividualSim(MapFormat,m,L,WholeArray,reactionprob)
	NormalisedData = Normalisation(AllPowerFirstSpecies,L,trials)
	return NormalisedData
end
