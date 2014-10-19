function output = readTimecourses( files )
% Runs readTimecourse over a group of files (separately).
% see readTimecourse.m
%

if nargin<1 || isempty(files),
    files = getFiles('*titration.txt');
end

output = [];

for i=1:numel(files)
    output(:,i) = readTimecourse( files{i} ); 
end

