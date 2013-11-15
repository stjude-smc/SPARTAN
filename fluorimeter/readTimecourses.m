function readTimecourses( files )
% Runs readTimecourse over a group of files (separately).
% see readTimecourse.m
%

if nargin<1 || isempty(files),
    files = getFiles('*titration.txt');
end

for i=1:numel(files)
    readTimecourse( files{i} ); 
end

