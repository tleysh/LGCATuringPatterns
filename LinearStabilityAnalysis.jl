#Stability Analysis module
module LinearStabilityAnalysis
    using SymEngine
    using DataFrames
    using CSV
    using Plots


    #Function to reformat the map expression
    function Reformatting(digs::Vector{Int64})
        extros = 16 - length(digs)
        extras = convert(Array{Int64},zeros(extros))
        revdigs = reverse(digs)
        alltogthertwo = vcat(extras,revdigs)
        return alltogthertwo
    end


    #Indicator functions for the expression of velocity channel occupation
    function Indicators(x::String)
        vars = [SymEngine.symbols("$x"*"1"),SymEngine.symbols("$x"*"2"),SymEngine.symbols("$x"*"3")]
        Indicator = Array{SymEngine.Basic}(4);
        Indicator[1] = (1-vars[1])*(1-vars[2])*(1-vars[3])
        Indicator[2] = ((1-vars[1])*(1-vars[2])*vars[3]+(1-vars[2])*(1-vars[3])*vars[1]+(1-vars[3])*vars[2]*(1-vars[1]))
        Indicator[3] = ((1-vars[1])*vars[2]*vars[3]+(1-vars[2])*vars[1]*vars[3]+(1-vars[3])*vars[2]*vars[1])
        Indicator[4] = (vars[1]*vars[2]*vars[3])
        return Indicator
    end


    function SymbolicDifferenceExpression(XExpression::Vector{SymEngine.Basic},YExpression::Vector{SymEngine.Basic},FirstSpeciesVector::Vector{Int64},SecondSpeciesVector::Vector{Int64})
        FullSymbolicExpressionX = SymEngine.Basic(0) ; FullSymbolicExpressionY = SymEngine.Basic(0);
        for j = 1:4
            for i = 1:4
                BothSpeciesComponent = XExpression[j]*YExpression[i]; inde = (j-1)*4+i;
                FirstSpeciesCoefficeint = FirstSpeciesVector[inde];
                SecondSpeciesCoefficeint = SecondSpeciesVector[inde];
                acompone = FirstSpeciesCoefficeint*BothSpeciesComponent/3;
                acomptwo = SecondSpeciesCoefficeint*BothSpeciesComponent/3;
                FullSymbolicExpressionX = FullSymbolicExpressionX+acompone ; FullSymbolicExpressionY = FullSymbolicExpressionY+acomptwo ;
            end
        end
        FullSymbolicExpressionX = FullSymbolicExpressionX-SymEngine.symbols("x1"); FullSymbolicExpressionY = FullSymbolicExpressionY-SymEngine.symbols("y1");
        return FullSymbolicExpressionX,FullSymbolicExpressionY
    end


    function ProducingJacobianElements(VariablesX::Vector{SymEngine.Basic},VariablesY::Vector{SymEngine.Basic},x::String,SteadyStateForX::Float64,SteadyStateForY::Float64,DifferenceExpression::SymEngine.Basic)
        JacobianElement = [];
        push!(JacobianElement,expand(DifferenceExpression));
        for i = 1:length(VariablesX)
            JacobianElementer = SymEngine.subs(JacobianElement[i],VariablesX[i],SteadyStateForX)
            push!(JacobianElement,JacobianElementer)
        end
        adder = length(VariablesX)
        for i = 1:length(VariablesY)
            JacobianElementer = SymEngine.subs(JacobianElement[i+adder],VariablesY[i],SteadyStateForY)
            push!(JacobianElement,JacobianElementer)
        end
        thediffer = SymEngine.symbols("$x"*"1")
        JacobianElement3 = expand(JacobianElement[end])
        JacobianElement4 = diff(JacobianElement3,thediffer)
        return JacobianElement4
    end

    # Jacobian
    function NewWholeJac(DifferenceExpressionX::SymEngine.Basic,DifferenceExpressionY::SymEngine.Basic,steds::Array{Float64,1},x1::SymEngine.Basic,x2::SymEngine.Basic,x3::SymEngine.Basic,y1::SymEngine.Basic,y2::SymEngine.Basic,y3::SymEngine.Basic)
        dfdx = ProducingJacobianElements([x2,x3],[y1,y2,y3],"x",steds[1],steds[2],DifferenceExpressionX); dfdy = ProducingJacobianElements([x1,x2,x3],[y2,y3],"y",steds[1],steds[2],DifferenceExpressionX);
        dgdx = ProducingJacobianElements([x2,x3],[y1,y2,y3],"x",steds[1],steds[2],DifferenceExpressionY); dgdy = ProducingJacobianElements([x1,x2,x3],[y2,y3],"y",steds[1],steds[2],DifferenceExpressionY);
     return dfdx,dfdy,dgdx,dgdy
    end

    # Eigns absolute value
    function EignsAbsoluteValues(ma::Int64,mi::Int64,q::Int64,w1::Float64,w2::Float64,w3::Float64,w4::Float64)
        uA = (1 + 2*cos(2*π*ma*q/101)); uI = (1+ 2*cos(2*π*mi*q/101));
        λ1 = 0.5*(w1*uA+w4*uI+sqrt(4*(w2*w3-w1*w4)*uA*uI+(w1*uA+w4*uI)^2));
        λ2 = 0.5*(w1*uA+w4*uI-sqrt(4*(w2*w3-w1*w4)*uA*uI+(w1*uA+w4*uI)^2));
        return abs(λ1),abs(λ2)
    end
    #Eigns with absolute, real and imaginary parts
    function EignsAbsoluteRealImagValues(ma::Int64,mi::Int64,q::Int64,w1::Float64,w2::Float64,w3::Float64,w4::Float64)
        uA = (1 + 2*cos(2*π*ma*q/101)); uI = (1+ 2*cos(2*π*mi*q/101));
        λ1 = 0.5*(w1*uA+w4*uI+sqrt(4*(w2*w3-w1*w4)*uA*uI+(w1*uA+w4*uI)^2));
        λ2 = 0.5*(w1*uA+w4*uI-sqrt(4*(w2*w3-w1*w4)*uA*uI+(w1*uA+w4*uI)^2));
        #return λ1,λ2
        return abs(λ1),abs(λ2),real(λ1),imag(λ1),real(λ2),imag(λ2)
    end
