#testing parallel programming
#Making optimisation changes
module SteadyStateFindingCore
    using SymEngine
    using NLsolve
    using CSV
    using DataFrames

    #indicator functions will be symbolic
    function IndicatorFunction(x::String)
        x = SymEngine.symbols("$x")
        Indicator = Array{SymEngine.Basic}(4);
        Indicator[1] = ((1-x)*(1-x)*(1-x))
        Indicator[2] = (3*(1-x)*(1-x)*x)
        Indicator[3] = (3*(1-x)*x*x)
        Indicator[4] = (x*x*x)
        return Indicator
    end


    # Function to convert map into symbolic difference equation
    function SymbolicDifferenceExpression(IndicatorFuncX::Vector{SymEngine.Basic},IndicatorFuncY::Vector{SymEngine.Basic},StandardTuringone::Vector{Int64},StandardTuringtwo::Vector{Int64})
        togetherone = SymEngine.Basic(0) ; togethertwo = SymEngine.Basic(0);
        for j = 1:4
            for i = 1:4
                holdtogether = IndicatorFuncX[j]*IndicatorFuncY[i]; inde = (j-1)*4+i;
                factone = StandardTuringone[inde];
                facttwo = StandardTuringtwo[inde];
                acompone = factone*holdtogether/3;
                acomptwo = facttwo*holdtogether/3;
                togetherone = togetherone+acompone ; togethertwo = togethertwo+acomptwo ;
            end
        end
        togetherxx = togetherone-SymEngine.symbols("x[1]"); togetheryy = togethertwo-SymEngine.symbols("x[2]");

        return togetherxx,togetheryy
    end


    # function that does findin and spits out true or false
    function Finder(sted::Vector{Float64},steds::Vector{Array{Float64,1}})
        bull = false
        for i in steds
            if sted == i
                bull = true
                break
            end
        end
        return bull
    end


    #Function that locates the steady states numerically be using the nlsolve package and 5 starting postions between 0 and 1
    function  NumericalSteadyStateLocator(MetaStringSymbolicExpression!)
        steds = Vector{Float64}[]
        for i in [[0.5,0.5],[0.25,0.25],[0.75,0.75],[0.25,0.75],[0.75,0.25],[0.1,0.1],[0.9,0.9]]
            #produce steady states from random entires between 0 and 1.
            nlsolveResults = nlsolve(MetaStringSymbolicExpression!,i)
            # Check if the steady states have converged
            if converged(nlsolveResults) == true
                    #if true then find the rest of the 1000 that are the same
                xsteadystate = round(nlsolveResults.zero[1],5)
                ysteadystate = round(nlsolveResults.zero[2],5)
                answer = Finder([xsteadystate,ysteadystate],steds)
                    if answer == false
                        if  (0.01 < xsteadystate < 0.99) && (0.01 < ysteadystate < 0.99)
                        #check they are within the right range
                            push!(steds,[xsteadystate,ysteadystate]);
                        end
                    end
            end
            # The length of the x and y arrays tells us the number of different steady states
        end
        numofsteds = length(steds);
        return steds,numofsteds
    end


    #Function for providing meta-string of the system of difference equations
    function ProducingMetaStringSymbolicExpression(XStringExpression::String,YStringExpression::String, mapnum::Int64)
        eval(parse("function test$mapnum(F,x); F[1] ="*XStringExpression*"; F[2] = "*YStringExpression*"; end;"))
    end

    #Function for providing meta-string of the system of difference equations
    function ProducingMetaStringSymbolicExpressionSingle(XStringExpression::String,YStringExpression::String, mapnum::String)
        #this is the function used for finding the steady states
        eval(parse("function test$mapnum(F,x); F[1] ="*XStringExpression*"; F[2] = "*YStringExpression*"; end;"))
    end

    #Function to reformat the map expression
    function Reformatting(digs::Vector{Int64})
        extros = 16 - length(digs)
        extras = convert(Array{Int64},zeros(extros))
        revdigs = reverse(digs)
        alltogthertwo = vcat(extras,revdigs)
        return alltogthertwo
    end

    #Function for parsing a full topology into Symbolic String expressions
    function ProducingalltheMetaStringSymbolicExpressions(j::Int64,batchnum::Int64,mapsperbatch::Int64)
        IndicatorFuncX = IndicatorFunction("x[1]"); IndicatorFuncY = IndicatorFunction("x[2]");
        elementos = DataFrame(CSV.read("/Users/tom/Documents/TopologyMapDictionaries/MapnumberdictionaryforAllTop$j.csv"))[2:3]
        starto = batchnum*mapsperbatch+1
        endo = (batchnum+1)*mapsperbatch
        @time for i = starto:endo
            print("Doing map $i \n")
            fullone = Reformatting(digits(elementos[1][i]))
            fulltwo = Reformatting(digits(elementos[2][i]))
            SymbolicX,SymbolicY = SymbolicDifferenceExpression(IndicatorFuncX,IndicatorFuncY,fullone,fulltwo)
            SymbolicStringX = string(SymbolicX); SymbolicStringY = string(SymbolicY);
            ProducingMetaStringSymbolicExpression(SymbolicStringX,SymbolicStringY,i)
        end
        return
    end

    function ProducingOneMetaStringSymbolicExpression(fullone,fulltwo,nameo)
        IndicatorFuncX = IndicatorFunction("x[1]"); IndicatorFuncY = IndicatorFunction("x[2]");
        SymbolicX,SymbolicY = SymbolicDifferenceExpression(IndicatorFuncX,IndicatorFuncY,Reformatting(digits(fullone)),Reformatting(digits(fulltwo)))
        SymbolicStringX = string(SymbolicX); SymbolicStringY = string(SymbolicY);
        ProducingMetaStringSymbolicExpressionSingle(SymbolicStringX,SymbolicStringY,nameo)
        return
    end

    #Parsing and evaluating one map
    function OneMap(mapnum::Int64)
        steds,numofdiffSteds = NumericalSteadyStateLocator(eval(parse("test$mapnum")))

        return steds,numofdiffSteds
    end

    #For topologies with large numbers of maps split into optimum batches of 1000 to not overflow the workspace
    #Use a shell script to go through full topology 
    function OneBatch(batchnum::Int64,Toponumber::Int64,mapsperbatch::Int64)
        df = DataFrame(mapID = Int64[], steds = String[])
        starto = batchnum*mapsperbatch+1
        endo = (batchnum+1)*mapsperbatch
        for i = starto:endo
            steds,numofdiffSteds = OneMap(i)
            if numofdiffSteds != 0
                push!(df,(i,string(steds)[17:end]))
                #push!(df,(i,string(steds)))
            end
            print("Doing map $i \n ")
        end
        CSV.write("/Users/tom/Documents/AsymptoticOctober/Topology$Toponumber"*"MondayBatchTest$batchnum.csv",df)
        return
    end
end

# -------------
#ITERATE THROUGH BASH NUMBER -
#Toponum = 15 ; Batchnum = 0 ; mapsperbatch = 1000 ; #default 1000 per a batch   941192
#Toponumber = Toponum+1;
#Toponumber = 15 ;
#ProducingalltheMetaStringSymbolicExpressions(Toponumber,Batchnum,mapsperbatch)
#print("\n This is the toponumber : $Toponumber \n")
#print("\n This is the batchnum : $Batchnum \n")
#print("\n This is the mapsperbatch : $mapsperbatch \n")
#@time OneBatch(Batchnum,Toponumber,mapsperbatch)
