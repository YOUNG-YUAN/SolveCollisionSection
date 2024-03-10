module PotentialFunction
    using   ..Units
    export  SetPotential, PFunc

    function SetPotential(FunctionType::Symbol, a=0.0::Float64, b=0.0::Float64, c=0.0::Float64)  # Set the function of Potential

        global PFunc  # Create closure (创建闭包函数，这使得它返回的是一个函数，而不是一个值)
        # You can use PFunc to Calculate the value of potential by specified, SetPotential should be done before this.
        if      FunctionType == :LennardJonesPotential  # Lennard-Jones Potential
            PFunc = x -> a * ((b / x).^12 - 2 * (b / x).^6)

        elseif FunctionType == :Free  # V = 0, free particle
            PFunc = x -> 0 * x

        elseif FunctionType == :FiniteDeepPotentialWell  # Finite deep potential well
            PFunc = x -> x >= a && x <= b ? c : 0

        elseif FunctionType == :HarmonicPotential  # Harmonic potential
            PFunc = x -> a * (x - b).^2

            # Other functions
        else 
            error("Potential function name is wrong! There is no such potential function.")
        end
        
    end

end