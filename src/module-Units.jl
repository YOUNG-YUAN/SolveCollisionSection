module Units

    export SetUnits, eV, meV, keV, nm, Å, Light, Me, ħ, ℎ, AddUnits

    function SetUnits(UnitsType::Symbol)
        
        if UnitsType == :SI     # SI 国际单位制
            # Warning, if numeric value is too small, which may cause a zero value, especially for exp() the function!
            # 注意，如果数值太小，可能导致零值，尤其在当你使用exp()时。
            global  eV     = 1.602176565e-19   # electron volt 电子伏特    [J]
            global  meV    = 1.602176565e-22
            global  keV    = 1.602176565e-16
            global  nm     = 1e-9              # nanometer 纳米            [m]
            global  Å      = 1e-10             # Angstrom 埃米             [m]
            global  Light  = 299792458         # speed of light            [m/s]
            global  Me     = 9.10938215e-31    # Mass of electron          [kg]
            global  ħ      = 1.05457266e-34    # reduced Plank constant    [J·s]
            global  ℎ      = ħ * 2 * pi        # Plank constant            [J·s]     
            
        elseif UnitsType == :Atom # 原子尺寸的单位制（并非原子单位制）
            global  eV      = 1.0                       # electron volt 电子伏特    [eV]
            global  meV     = 1.0e-3
            global  keV     = 1.0e+3
            global  nm      = 1.0e+1                    # nanometer 纳米            [Å]
            global  Å       = 1.0                       # Angstrom 埃米             [Å]
            global  Light   = 2.99792458e18             # speed of light            [Å/s]
            global  Me      = 5.10998885e+5*eV/Light^2  # Mass of electron          [Å/s]
            global  ħ       = 6.582119514e-16           # reduced Plank constant    [eV·s]
            global  ℎ       = ħ * 2 * pi                # Plank constant            [eV·s]
            
            # Customize Units (in future)
        end
    end



    function AddUnits(NewUnit::String)  # Will improve in future
        """
        Only for ACSII, UTF-8 is not supported!!! Will improve this function in future.
        """
        #=
        You can add your custom units by this function
        Format of NewUnit is like "ρ = 2nm", Double quotation marks are needed, and an equation is between them.
        =#
        if length(findall("=",str)) == 1  # Only one equal sign——"="
            eval(Meta.parse("global " * NewUnit))   # create NewUnit as global variables
            UnitName = NewUnit[1 : findnext("=", NewUnit, 1)[1] - 1]  # get the UnitName from NewUnit the string
            eval(Meta.parse("export " * UnitName))  # export the new unit
        else
            error("NewUnit format is incorrect!")
        end
    end
end
