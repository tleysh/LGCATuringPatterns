#Core simulator
module SimulationFunctions
    using CSV
    using DataFrames

    #Function for reformatting the map inputs in two {Int8}Arrays
    function reformating(SpeciesInput::Int64)
        firstdigs = digits(SpeciesInput)
        diggs = map(v->Int8(v),firstdigs)
        extros = 16 - length(diggs)
        extras = map(v -> Int8(v),zeros(extros))
        vcat(extras,reverse(diggs))
        revdigs = reverse(diggs)
        alltogthertwo = vcat(extras,revdigs)
        return alltogthertwo
    end

    #Function for making the reformatted arrays into a combined map
    function makingMap(newfullone::Array{Int8,1},newfulltwo::Array{Int8,1})
             themappa = Array{Int8,2}(undef,16,2)
             for i = Int8(1):Int8(16)
                  themappa[i,1] = newfullone[i]
                  themappa[i,2] = newfulltwo[i]
             end
             return themappa
    end

    #Function to initialise the spatial-temporal domain
    function InitialArray(M::Int64,L::Int64)
        #M = this is the number of species
        # L is the length of array
        AllSpeciesone = rand([0,1],101,3)
        AllSpeciestwo = rand([0,1],101,3)
        AllSpeciesone8bit = map(v->Int8(v),AllSpeciesone)
        AllSpeciestwo8bit = map(v->Int8(v),AllSpeciestwo)
        AllSpecies = [AllSpeciesone8bit,AllSpeciestwo8bit]
        return AllSpecies
    end

    #----- Reaction functionality -----

    #Function to Combine the individual velocity channels ready for the reaction step
    function Combining8Bitchannels(IndividualChannelDomain::Array{Array{Int8,2},1},L::Int64)
        combinedvelocitychannelsFirst = map(v->Int8(v),[sum(IndividualChannelDomain[1][x,1:3]) for x in 1:L])
        combinedvelocitychannelsSecond = map(v->Int8(v),[sum(IndividualChannelDomain[2][x,1:3]) for x in 1:L])
        return [combinedvelocitychannelsFirst combinedvelocitychannelsSecond]
    end

    #Reaction step returns Combined Channel format
    function ReactionStep(Mappo::Array{Int8,2},L::Int64,CombinedChannels::Array{Int8,2},probabilitythresholdfirstspecies::Float64,probabilitythresholdsecondspecies::Float64)
        Allp = rand(L); Allq = rand(L);
        afterreaction = deepcopy(CombinedChannels)
        for x in 1:L
            if Allp[x] <= probabilitythresholdfirstspecies
                n,m = CombinedChannels[x,1:2]
                pos = n*4+(m+1)
                updatedn,updatedm = Mappo[pos,1:2]
                afterreaction[x,1] = updatedn
            end
            if Allq[x] <= probabilitythresholdsecondspecies
                n,m = CombinedChannels[x,1:2]
                pos = n*4+(m+1)
                updatedn,updatedm = Mappo[pos,1:2]
                afterreaction[x,2] = updatedm
            end
        end
        return afterreaction
    end


    #----- Diffusion functionality -----

    #The possible states the 3 channel model can be within
    function ThePossibleChannelArrangements()
        Array1 = [[Int8(0) Int8(0) Int8(0)]]; Array2 = [[Int8(1) Int8(0) Int8(0)],[Int8(0) Int8(1) Int8(0)],[Int8(0) Int8(0) Int8(1)]]; Array3 = [[Int8(1) Int8(1) Int8(0)],[Int8(0) Int8(1) Int8(1)],[Int8(1) Int8(0) Int8(1)]]; Array4 = [[Int8(1) Int8(1) Int8(1)]];
        WholeArray = [Array1,Array2,Array3,Array4];
        return WholeArray
    end

    #8 bit splitter works on a single species and chooses a new arrangement of the channels for combined channels of a single species.
    function splitter8bit(Aspecies::Int8,WholeArray::Array{Array{Array{Int8,2},1},1})

        PossibleChoices = WholeArray[Aspecies+1]
        if (Aspecies+1) == 1 ||  (Aspecies+1) == 4
            theArrangement = 1;
        else
            theArrangement = rand([1,2,3])
        end
        return PossibleChoices[theArrangement]
    end


    #SplitTogether take two species and splits these up by apply splitter8bit to each species seperately
    function SplitTogether(TwoSpecies::Array{Int8,1},WholeArray::Array{Array{Array{Int8,2},1},1})
            splitfirst,splitsecond = splitter8bit(TwoSpecies[1],WholeArray),splitter8bit(TwoSpecies[2],WholeArray)
            return [splitfirst,splitsecond]
    end


    #AllSplitTogether takes the two species and returns the individual channel domains
    function AllSplitTogether(CombinedChannelsDomain::Array{Int8,2},WholeArray::Array{Array{Array{Int8,2},1},1},L::Int64)
            TheSplitnewone = Array{Int8,2}(undef,101,3)
            TheSplitnewtwo = Array{Int8,2}(undef,101,3)
            for x in 1:L
                TheSplitnewone[x,1:3],TheSplitnewtwo[x,1:3] = SplitTogether(CombinedChannelsDomain[x,1:2],WholeArray)
            end
            return [TheSplitnewone,TheSplitnewtwo]
    end


    #Diffusion Jumper moves the channel values by the predetermined diffusion rate.
    function Jumper8bitAttempttwo(AfterSplit,m::Array{Int64,1},L::Int64)
            return [[vcat(AfterSplit[1][m[1]+1:L,1],AfterSplit[1][1:m[1],1]) AfterSplit[1][1:L,2] vcat(AfterSplit[1][L-m[1]+1:L,3],AfterSplit[1][1:L-m[1],3])],[vcat(AfterSplit[2][m[2]+1:L,1],AfterSplit[2][1:m[2],1]) AfterSplit[2][1:L,2] vcat(AfterSplit[2][L-m[2]+1:L,3],AfterSplit[2][1:L-m[2],3])]]
    end

    # Diffusion Step contains both the splitting (shuffling) step and also the jump diffusion
    function DiffusionStep8bit(afterreaction::Array{Int8,2},m::Array{Int64,1},L::Int64,WholeArray::Array{Array{Array{Int8,2},1},1})
            #Split up into channels after reaction
            AfterSplit =  AllSplitTogether(afterreaction,WholeArray,L)
            #Jump channels across the lattice
            #afterdiffusion = Jumper8bit(AfterSplit,m,L)
            AfterDiffusion = Jumper8bitAttempttwo(AfterSplit,m,L)
            return AfterDiffusion
    end

    # ---- One Single Time-Step -------

    #One step involves a reaction step and then a diffusion step
    function Onestep8bit(IndividualChannelDomain::Array{Array{Int8,2},1},m::Array{Int64,1},Mappo::Array{Int8,2},L::Int64,probabilitythresholdfirstspecies::Float64,probabilitythresholdsecondspecies::Float64,WholeArray::Array{Array{Array{Int8,2},1},1})
        #Reactionstep
        CombinedChannelsDomain = Combining8Bitchannels(IndividualChannelDomain,L)
        AfterReaction = ReactionStep(Mappo,L,CombinedChannelsDomain,probabilitythresholdfirstspecies,probabilitythresholdsecondspecies)
        AfterDiffusion = DiffusionStep8bit(AfterReaction,m,L,WholeArray)

        return AfterDiffusion
    end
end
