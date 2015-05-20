function[data, DETERMINE_strokeVolumes, DETERMINE_ejectionFractions,  MESA_strokeVolumes,  MESA_ejectionFractions] = calcEjectionFraction(data)

for i = 1:401
    %   i
    data(i).strokeVolume = data(i).diastolic.endo.volume - data(i).systolic.endo.volume;
    data(i).ejectionFraction = ((data(i).strokeVolume)/data(i).diastolic.endo.volume)*100;
end

for n = 1:100
for d = data(1).DETERMINE_indices(n)
    
    data(d).DETERMINE.strokeVolume =  data(d).strokeVolume;
    data(d).DETERMINE.ejectionFraction = data(d).ejectionFraction;
    
    DETERMINE_strokeVolumes(n) =  data(d).strokeVolume;
    DETERMINE_ejectionFractions(n) = data(d).ejectionFraction;
end

for m = data(1).MESA_indices(n)
    
    data(m).MESA.strokeVolume = data(m).strokeVolume;
    data(m).MESA.ejectionFraction = data(m).ejectionFraction;
    
    MESA_strokeVolumes(n) = data(m).strokeVolume;
    MESA_ejectionFractions(n) = data(m).ejectionFraction;
    
end

end
end
