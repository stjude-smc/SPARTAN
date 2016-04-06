function output = settingsDialog(input,varargin)
%function output = settingsDialog(input, fields,prompt,types, fun,fin)
%settingsDialog  Input dialog to change arbitrary struct fields
%
%   OPT = settingsDialog(OPT) displays an inputdlg with one line for each
%   field in the struct OPT, allowing a user to easily change the values.
%
%   OPT = settingsDialog(OPT,FIELDS,PROMPTS,TYPES) specifies which FIELDS in 
%   OPT to show, custom dialog PROMPTS, and TYPES for input validation.
%   Each of these cell arrays has the same size. TYPES are class() types like
%   'double' or can be a cell array of allowed values.
%
%   settingsDialog(...,FUN,FIN) after the dialog closes, calls the function
%   handle as FUN(OPT,FIN{:}).

%   Copyright 2007-2016 Cornell University All Rights Reserved.



%% Process input arguments
nargoutchk(0,1);
narginchk(1,6);
if ~isstruct(input) || numel(input)>1,
    error('First argument must be a scalar struct');
end

% Parse out field and function arguments from beginning and end of varargin.
idxFun = find(  cellfun( @(x)isa(x,'function_handle'), varargin )  );
assert( numel(idxFun)<=1, 'Multiple function inputs now allowed' );
if ~isempty(idxFun)
    args = varargin(1:idxFun-1);
else
    args = varargin;
end
    

% Construct default field names, prompts, etc.
fields = fieldnames(input);
prompt = strcat(fields,':');
types = cell(numel(fields),1);

switch numel(args)
    case 0
    case 1
        fields = args{1};
    case 2
        [fields,prompt] = args{:};
    case 3
        [fields,prompt,types] = args{:};
    otherwise
        error('Too many input arguments');
end


%% Prompt for new values
output = input;

% Select parameters to show, converting numbers to strings for inputdlg().
currentopt = cellfun( @(x)num2str(input.(x)), fields, 'UniformOutput',false );

% Prompt user for new values
answer = inputdlg(prompt, [mfilename ' display settings'], 1, currentopt);
if isempty(answer), return; end  %user hit cancel

% Reformat user answers to original format.
for k=1:numel(answer),
    value = input.(fields{k});
    
    if isnumeric(value)
        value = str2double(answer{k});
        
    elseif islogical(value)
        switch lower(answer{k})
            case {'on','true','yes'}
                value = true;
            case {'off','false','no'}
                value = false;
            otherwise
                value = logical(str2double(answer{k}));
        end
        
    else
        value = answer{k};
    end
    
    % Enforce type constraints
    if iscell(types{k}) && ~ismember(value,types{k}),
        errordlg( { sprintf('Invalid value for %s. Must be one of:',fields{k}); ...
                    strjoin(types{k},', ') } );
        output = input;
        return;
    end
    
    output.(fields{k}) = value;
end



%% Execute callback to update target with new settings.
if ~isempty(idxFun),
    fun  = varargin{idxFun};
    fin  = varargin{idxFun+1};
    fun(fin{:}, output);
end


end %function settingsDialog


