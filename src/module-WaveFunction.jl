module WaveFunction 

    using ..PotentialFunction, ..Units

    export funcpoint, WaveFunctionSetting

    export SetWaveFunction, SolveWaveFunction, globalWFSetting, ShootingMethod

    mutable struct funcpoint
        R           ::Float64
        Value       ::Float64
    end

    mutable struct WaveFunctionSetting
        M                   ::Float64
        Rmin                ::Float64
        Rmax                ::Float64
        Step                ::Float64
        iterativeDirection  ::Symbol
        initialValue1       ::Float64
        initialValue2       ::Float64
    end

    function SetWaveFunction(M:: Float64, Rmin::Float64, Rmax::Float64, Step::Float64, iterativeDirection::Symbol, initialValue1::Float64, initialValue2::Float64)
        # Set the range, iterative step, iterative initial two points, iterative direction of wavefunction.
        global globalWFSetting
        globalWFSetting = WaveFunctionSetting(M, Rmin, Rmax, Step, iterativeDirection, initialValue1, initialValue2)
    end

    function SolveWaveFunction(L::Int64, Energy::Float64, WaveFunctionSettings::WaveFunctionSetting=globalWFSetting)
        # if !(globalWFSetting in varnames(Main))
        #     error("globalWFSetting does not exist! Please use SetWaveFunction() before using SolveWaveFunction()!")
        # end
        M                   = WaveFunctionSettings.M
        Rmin                = WaveFunctionSettings.Rmin
        Rmax                = WaveFunctionSettings.Rmax
        Step                = WaveFunctionSettings.Step
        iterativeDirection  = WaveFunctionSettings.iterativeDirection
        initialValue1       = WaveFunctionSettings.initialValue1
        initialValue2       = WaveFunctionSettings.initialValue2
        
        # Set the initial condition of wavefunction, at least two points.
        """
        Numerov method
        The radial Schrodinger equation is :
                Eu(r) = -ħ²/2m * d²u(r)/dr² + L(L+1)ħ²/2mr² * u(r) + V(r)u(r) .
        Now we do some change, the radial Schrodinger equation can read:
                -ħ²/2m * d²u(r)/dr² = F(L,r,E)u(r) ,
                F(L,r,E) = V(r) + L(L+1)ħ²/2mr² .
        According to Numerov method, radial wavefunction can solve by:
                w(r + step) = 2w(r) - w(r - step) + step^2 * F(L,r,E)u(r)
                u(r) = [1 - step^2 / 12 * F(L,r,E)]^(-1) * w(r)
        I call [1 - step^2 / 12 * F(L,r,E)] as exchange coefficient
        We just need two initial points, then we can solve radial wavefunction
        """

        function F_LrE(M::Float64, L::Int64, r::Float64, Energy::Float64)  
            F = L * (L + 1) ./ r.^2 + 2 * M / ħ^2 * (- Energy + PFunc(r))
            return F
        end

        function ExchangeCoeff(Step::Float64, r::Float64, L::Int64, Energy::Float64, M::Float64)  # It derives from Numerov method
            if r > 0.0  # r > 0 is valid
                ec = 1 - Step^2 * F_LrE(M, L, r, Energy) / 12
            elseif r == 0.0  # r = 0 is a singular point
                ec = 0
            else  # r < 0 is invalid
                error("ERROR: Wrong r of RadialWaveFunc!")
            end
            return ec
        end

        acomRWF = Vector{funcpoint}()  # acompany Radial WaveFunction w(r)
        RWF = Vector{funcpoint}()  # Radial WaveFunction u(r)
        WF = Vector{Float64}()  # WaveFunction ϕ(r), is for normalization

        if iterativeDirection == :right  # iteration is from the left side to the right side

            for i = Int64(ceil(Rmin/Step)) : 1 : Int64(ceil(Rmax/Step))
                # Caution: i = 0 corresponds to R = 0, but index = 1 for acomRWF the vector !
                j = i + 1 - Int64(ceil(Rmin/Step))  # j is the index of acomRWF where acomRWF.R[j] = i * step
                if j >= 3
                    value       = 2 * acomRWF[j-1].Value - acomRWF[j-2].Value + Step^2 * 
                                    F_LrE(M,L,(i-1)*Step,Energy) * RWF[j-1].Value
                    RWFvalue    = value / ExchangeCoeff(Step, i*Step, L, Energy, M)
                    RWFpoint    = funcpoint(i * Step, RWFvalue)
                    push!(acomRWF, funcpoint(i * Step, value))
                    push!(RWF, funcpoint(i * Step, RWFvalue))
                    # push!(WF, RWFvalue / (i*Step))
                    
                elseif j == 2
                    value       = initialValue2 * ExchangeCoeff(Step, i*Step, L, Energy,M)
                    aRWFpoint   = funcpoint(i * Step, value)
                    RWFpoint    = funcpoint(i * Step, initialValue2)
                    push!(acomRWF, funcpoint(i * Step, value))
                    push!(RWF, funcpoint(i * Step, initialValue2))
                    # push!(WF, initialValue2 / (i*Step))

                elseif j == 1
                    value       = initialValue1 * ExchangeCoeff(Step, i*Step, L, Energy, M)
                    aRWFpoint   = funcpoint(i * Step, value)
                    RWFpoint    = funcpoint(i * Step, initialValue1)
                    push!(acomRWF, funcpoint(i * Step, value))
                    push!(RWF, funcpoint(i * Step, initialValue1))
                    # push!(WF, initialValue1 / (i*Step))
                end
            end

        elseif iterativeDirection == :left  # iteration is from the right side to the left side
            for i = Int64(ceil(Rmax/Step)) : -1 : Int64(ceil(Rmin/Step))
                # Caution: i = 0 corresponds to R = 0, but index = 1 for acomRWF the vector !
                j = Int64(ceil(Rmax/Step)) - i + 1  # j is the index of acomRWF where acomRWF.R[j] = i * step
                if j >= 3
                    value       = 2 * acomRWF[j+1].Value - acomRWF[j+2].Value + Step^2 * 
                                    F_LrE(M,L,(i+1)*Step,Energy) * RWF[j+1].Value
                    aRWFpoint   = funcpoint(i * Step, value)
                    aRWFpoint   = funcpoint(i * Step, value)
                    RWFvalue    = value / ExchangeCoeff(Step, i*Step, L, Energy, M)
                    RWFpoint    = funcpoint(i * Step, RWFvalue)
                    push!(acomRWF, funcpoint(i * Step, value))
                    push!(RWF, funcpoint(i * Step, RWFvalue))
                    # push!(WF, RWFvalue / (i*Step))

                elseif j == 2
                    value       = initialValue2 * ExchangeCoeff(Step, i*Step, L, Energy, M)
                    aRWFpoint   = funcpoint(i * Step, value)
                    RWFpoint    = funcpoint(i * Step, initialValue2)
                    push!(acomRWF, funcpoint(i * Step, value))
                    push!(RWF, funcpoint(i * Step, initialValue2))
                    # push!(WF, initialValue2 / (i*Step))

                elseif j == 1
                    value       = initialValue1 * ExchangeCoeff(Step, i*Step, L, Energy, M)
                    aRWFpoint   = funcpoint(i * Step, value)
                    RWFpoint    = funcpoint(i * Step, initialValue1)
                    push!(acomRWF, funcpoint(i * Step, value))
                    push!(RWF, funcpoint(i * Step, initialValue1))
                    # push!(WF, initialValue1 / (i*Step))
                end
            end

        else
            error("Wrong iterativeDirection!!")
        end

        # NormalizationCoeff = 1 / sum(abs2.(WF))
        # for i = eachindex(RWF)
        #     RWF[i].Value = NormalizationCoeff * RWF[i].Value  # Normalization
        # end

        return RWF
    end

    # Other functions for this module

end
