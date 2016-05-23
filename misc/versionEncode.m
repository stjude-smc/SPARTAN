function output = versionEncode(string)
% Convert version number string (x.y.z) to an integer for easy comparison.
% 2.12.6 becomes 2012006.

if ischar(string)
    version = cellfun( @(s)sscanf(s,'%f'), strsplit(string,'.') );
elseif isnumeric(string)
    version = string;
else
    error('Invalid input');
end

if numel(version)<3,
    version = [version zeros(1,max(0,3-numel(version)))];
elseif numel(version)>3
    error( ['Invalid version number ' string] );
end

output = 1e6*version(1) + 1e3*version(2) + version(3);

end
