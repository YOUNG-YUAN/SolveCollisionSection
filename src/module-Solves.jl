module Solves

    #include("module-WaveFunction.jl")

    using   ..Units, ..PotentialFunction, ..WaveFunction
    #using .PotentialFunction, .WaveFunction
    export GetPhaseShift, GetGrossCrossSection

    # mutable struct funcpoint
    #     x::Float64
    #     y::Float64
    # end

    function GetPhaseShift(RWF::Vector{funcpoint}, Rv::Float64, L::Int64, Energy::Float64, WaveFunctionSettings::WaveFunctionSetting=globalWFSetting)

        M       = WaveFunctionSettings.M
        Step    = WaveFunctionSettings.Step

        # Calculate the phase of collision
        R = Vector{Float64}()
        Wave = Vector{Float64}()
        for i in eachindex(RWF)
            append!(R, RWF[i].R)
            append!(Wave, RWF[i].Value)
        end
        """
        When r > rmax, where rmax is the margin of V(r), radial wavefunction ul is proportion to following formula:
                ul(r > rmax) ∝ kr[jl(kr) * cos(δl) - nl(kr) * sin(δl)] ,
                k = √(2mE/ħ²) ,
        where jl and nl are regular and irregular spherical Bessel functions respectively. If we find two point 
        beyond rmax, r1 and r2, then we can solve phase shift by following formula:
                tan(δl) = (Kjl(r1) - jl(r2)) / (Knl(r1) - nl(r2)) ,
                K = r1*ul(r2) / (r2*ul(r1)) .
        We can let r1 equals to rmax, but we need around it up to the grid. r2 is larger than r1 and it is aadvisable
        to take it roughly half a de Broglie wavelength beyond the latter, λ = 2π/k = 2πħ/√(2mE).
        For Lennard-Jones potaintial, a good value is rmax ≈ 5ρ. 

        """
        k = sqrt(2 * M * Energy / ħ^2)
        R1 = Step * ceil(Rv / Step)  # 5 means max margin of V(r) is 5ρ, we let r1 equals to rmax, and around it up to the grid
        R2 = Step * floor((R1 + 0.5 * 2 * pi / k)/Step)  # r1
        u1 = Wave[argmin(abs.(R .- R1))]
        u2 = Wave[argmin(abs.(R .- R2))]
        # argmin() is for find the index of element of R, the element equals to R1, then we can                                                             
        K = R1 * u2 / (R2 * u1)
        jl1 = SphericalBesselFunction(L, k * R1, :regular)
        jl2 = SphericalBesselFunction(L, k * R2, :regular)
        nl1 = SphericalBesselFunction(L, k * R1, :irregular)
        nl2 = SphericalBesselFunction(L, k * R2, :irregular)
        PhaseShift = atan((K * jl1 - jl2) / (K * nl1 - nl2))
        return PhaseShift
    end

    function SphericalBesselFunction(L::Int64, x::Float64, option::Symbol=:regular)
        
        if option == :regular  # regular spherical bessel function
            SBF0 = sin.(x) ./ x  # SBF0 means the zeroth order of spherical bessel function (regular), then SBF1 means the first order, etc. Same way in irregular SBF
            SBF1 = sin.(x) ./ x.^2 - cos.(x)./x
        elseif option == :irregular  # irregular spherical bessel function
            SBF0 = - cos.(x) ./ x
            SBF1 = - cos.(x) ./ x.^2 - sin.(x)./x
        else
            error("ERROR: Wrong option of SphericalBesselFunction(L::Int64, x::Float64, option::Symbol:regular)!")
        end

        if L == 0
            return SBF0
        elseif  L == 1
            return SBF1
        elseif L >= 2
            SBF_Vector = zeros(L+1)
            SBF_Vector[1] = SBF0
            SBF_Vector[2] = SBF1
            for i = 3 : L+1
                SBF_Vector[i] = (2i - 3) ./ x * SBF_Vector[i-1] - SBF_Vector[i-2]
            end
        else
            error("Error: Wrong L!")
        end
        return SBF_Vector[L+1]
    end

    function GetGrossCrossSection(Rv::Float64, L::Vector{Int64}, Energy::Vector{Float64}, M::Float64)
        
        GCS = zeros(length(Energy))  # GCS is a vector of all gross cross section
        PS = zeros(length(Energy))  # PS is a vector of all phaseshift
        for i = eachindex(Energy)
            k = sqrt(2 * M * Energy[i] / ħ^2)
            for j = eachindex(L)
                RWF = SolveWaveFunction(L[j], Energy[i])
                PhaseShift = GetPhaseShift(RWF, Rv, L[j], Energy[i])
                GCS[i] = GCS[i] + 4pi / k^2 * (2L[j] + 1) * (sin(PhaseShift))^2
                PS[i] = PS[i] + PhaseShift
            end
        end
        return (GCS,PS)
    end


end