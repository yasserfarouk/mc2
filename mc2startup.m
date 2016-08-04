% This file is a part of the MC2 toolbox developed by Y. Mohammand and T. Nishida.
%Please do not remove this comment
%
% Using this file is governed by the license of MC2 which you can find in LICENSE.md
% 
% You can find more information about this toolbox here:
% - Yasser Mohammad and Toyoaki Nishida, "MC2: An Integrated Toolbox for Change, Causality, 
%   and Motif Discovery", 29th International Conference on Industrial, Engineering & 
%   Other Applications of Applied Intelligent Systems (IEA/AIE) 2016, pp. 128 -- 141.
% - Yasser Mohammad and Toyoaki Nishida, "Data Mining for Social Robotics", Springer 2016.
%

disp('MC^2 Startup: (c) Yasesr Mohammad 2009-2019 http://www.ii.ist.i.kyoto-u.ac.jp/~yasser')
mc2path = fileparts( mfilename('fullpath') );

mainpath = fullfile(mc2path, 'algorithms');
if exist(mainpath,'dir')
    addpath(genpath(mainpath));    
end

mainpath = fullfile(mc2path, 'learning');
if exist(mainpath,'dir')
    addpath(genpath(mainpath));    
end

mainpath = fullfile(mc2path, 'preprocessing');
if exist(mainpath,'dir')
    addpath(genpath(mainpath));    
end

mainpath = fullfile(mc2path, 'preprocessing/transform');
if exist(mainpath,'dir')
    addpath(genpath(mainpath));    
end

mainpath = fullfile(mc2path, 'helpers');
if exist(mainpath,'dir')
    addpath(genpath(mainpath)); 
end

mainpath = fullfile(mc2path, 'demos');
if exist(mainpath,'dir')
    addpath(genpath(mainpath));    
end


mainpath = fullfile(mc2path, 'generators');
if exist(mainpath,'dir')
    addpath(genpath(mainpath)); 
end

mainpath = fullfile(mc2path, 'distances');
if exist(mainpath,'dir')
    addpath(genpath(mainpath)); 
end

mainpath = fullfile(mc2path, 'eval');
if exist(mainpath,'dir')
    addpath(genpath(mainpath)); 
end

% mainpath = fullfile(mc2path, 'test');
% if exist(mainpath,'dir')
%     addpath(mainpath); 
% end

global MC2BINPATH;
MC2BINPATH=fullfile(mc2path, 'bin');

clear mainpath mc2path
disp('MC^2 Startup Completed.')