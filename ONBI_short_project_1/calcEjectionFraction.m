function[data] = calcEjectionFraction(data)

for i = 1:401
%   i
data(i).strokeVolume = data(i).diastolic.endo.volume - data(i).systolic.endo.volume;
data(i).ejectionFraction = (data(i).strokeVolume./data(i).diastolic.endo.volume)*100;
    
    
%     for d = 1:100
%         if data(1).DETERMINE_indices(d)==i
%             data(i).DETERMINE.strokeVolume = data(i).DETERMINE.diastolic.endo.volume - data(i).DETERMINE.systolic.endo.volume;
%             data(i).DETERMINE.ejectionFraction = (data(i).DETERMINE.strokeVolume./data(i).DETERMINE.diastolic.endo.volume)*100;
%         else
% %         end
%     end

%     for d = 1:100
%         if data(1).MESA_indices(d)==i
%             data(i).MESA.strokeVolume = data(i).MESA.diastolic.endo.volume - data(i).MESA.systolic.endo.volume;
%             data(i).MESA.ejectionFraction = (data(i).MESA.strokeVolume./data(i).MESA.diastolic.endo.volume)*100;
%         else
%         end
%     end
end