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