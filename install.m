% This files installs the MATLAB support package for AdaptiveStrucrure Design Methodology (AdaptiveStructres package).
% external dependencies calfem fem package 
% Copyright 2012 Gennaro Senatore
clear all;clc;close all; 
userpath('clear');
% look for main file
mainFolder=which('example.m','-all');

% make sure we are in the right folder and there are no other main.m files
if length(mainFolder) < 1, 
    msg=' Cannot find DR.m, please run this file from the folder containing main.m';
    error(msg);
end
if length(mainFolder) > 1,
    msg=' There is at least another main.m file in the path, please delete any other versions before installing this one';
    error(msg);
end

% get the main folder
mainFolder=mainFolder{1};mainFolder=mainFolder(1:end-9);
cd(mainFolder);
% Add target directories and save the updated path
addpath(genpath(mainFolder));
disp(' Structural Analysis folder added to the path');

% result = savepath;
% if result==1
%     nl = char(10);
%     msg = [' Unable to save updated MATLAB path (<a href="http://www.mathworks.com/support/solutions/en/data/1-9574H9/index.html?solution=1-9574H9">why?</a>)' nl ...
%            ' On Windows, exit MATLAB, right-click on the MATLAB icon, select "Run as administrator", and re-run install_arduino.m' nl ...
%            ' On Linux, exit MATLAB, issue a command like this: sudo chmod 777 usr/local/matlab/R2011a/toolbox/local/pathdef.m' nl ...
%            ' (depending on where MATLAB is installed), and then re open MATLAB and re-run install_arduino.m' nl ...
%            ];
%     error(msg);
% else
%     disp(' Saved updated MATLAB path');
%     disp(' ');
% end

%cd ..
%mainFolder=pwd;
modelDataFolder=strcat(mainFolder,'models\chain\');
addpath(fullfile(modelDataFolder,''));
%mainFolder=strcat(mainFolder,'\DR\');
%cd(fullfile(mainFolder,''));
clear wa result nl msg ap
%load fileBreakPoints;
%dbstop(breakPoints);
