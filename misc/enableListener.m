function output = enableListener(obj,value)
% Enable or disable event listener object, across MATLAB versions.
% enableListener( listenerHandle, binaryValue );
% listenerHandle may be an array.

if nargin<2,
    % Get the current value
    if ischar(obj(1).Enabled),  %MATLAB 2014a
        output = get(obj,'Enabled');

    elseif islogical(obj(1).Enabled),  %MATLAB 2015a
        output = [obj.Enabled];
    end

else
    % Set to new value
    if ischar(obj(1).Enabled),  %MATLAB 2014a
        set(obj,'Enabled',onoff(value));

    elseif islogical(obj(1).Enabled),  %MATLAB 2015a
        [obj.Enabled] = deal(value);
    end
end

end %FUNCTION enableListener