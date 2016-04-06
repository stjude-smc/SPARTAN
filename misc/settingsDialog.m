function opt = settingsDialog(opt, fields,prompt, fun, varargin)
%settingsDialog  Input dialog to change arbitrary struct fields
%
%   OPT = settingsDialog(OPT) displays an inputdlg with one line for each
%   field in the struct OPT, allowing a user to easily change the values.
%
%   OPT = settingsDialog(OPT,FIELDS,PROMPTS) specifies which fields in OPT to
%   show and optionally specify custom prompt text for each.
%
%   settingsDialog(...,FUN,FIN) calls the function handle FUN with a variable 
%   number of input arguments FIN once new values are chosen by the user.

%   Copyright 2007-2016 Cornell University All Rights Reserved.


nargoutchk(0,1);
assert(numel(opt)==1, 'Struct arrays not allows');

% Use structure field names if not given.
if nargin<2
    fields = fieldnames(opt);
end
if nargin<3
    prompt = strcat(fields,':');
end


% Get current parameter values
currentopt = cellfun( @(x)num2str(opt.(x)), fields, 'UniformOutput',false );

% Ask user for new ones
answer = inputdlg(prompt, [mfilename ' display settings'], 1, currentopt);
if isempty(answer), return; end  %user hit cancel

% Reformat user answers to original format.
for k=1:numel(answer),
    original = opt.(fields{k});
    
    if isnumeric(original)
        opt.(fields{k}) = str2double(answer{k});
    elseif islogical(original)
        opt.(fields{k}) = logical(str2double(answer{k}));
    else
        opt.(fields{k}) = answer{k};
    end
    
    % FIXME: verify string fields.
end


% Update target
if nargin>3,
    fun(varargin{:}, opt);
end


end %function settingsDialog


