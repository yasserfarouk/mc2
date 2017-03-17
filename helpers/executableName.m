function [ fullName ] = executableName( baseName )
%Gets the executable name adjusted to the machine type
%   

if ispc()
    fullName = strcat(baseName, '.exe');
elseif ismac()
    fullName = strcat(baseName, '.mac');
elseif isunix()
    fullName = baseName;
end


end

