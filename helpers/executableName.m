function [ fullName ] = executableName( baseName )
%Gets the executable name adjusted to the machine type
%   

if ispc()
    fullName = strcat('./', strcat(baseName, '.exe'));
elseif ismac()
    fullName = strcat('./', strcat(baseName, '.mac'));
elseif isunix()
    fullName = baseName;
end


end

