function [stringer,FLAGSTRING] = AddStringer()
% This function reads the stringer input from the Excel file.
% If no stringer input is find stringer is empty and FLAGSTRING = 0.

global xlsFileName

StringerData = xlsread(xlsFileName,'Stringers');

if isempty(StringerData)
    FLAGSTRING = 0;
    stringer = [];
else
    FLAGSTRING = 1;
    stringer.yloc = StringerData(:,1);
    for i=1:length(stringer.yloc)-1
        if stringer.yloc(i+1)==stringer.yloc(i) % Check whether two entries have the same y-location and offset by 1e-6 span to prevent errors in interpolation
            stringer.yloc(i+1) = stringer.yloc(i+1)+stringer.yloc(end)*1e-6;
        end
    end
    stringer.pitch = StringerData(:,3);
    stringer.EA = StringerData(:,4);
    stringer.mA = StringerData(:,5);
    stringer.height = StringerData(:,2);
end

end