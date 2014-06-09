classdef Traces < handle & matlab.mixin.Copyable
% Traces: abstract class defining the interface for Traces classes
%
%    You can create instances of this abstract class, but they will have no
%    data channels. Use TracesFret or TracesFret4, TracesFluor, etc. These have
%    a defined set of data channels, but use the methods in this superclass for
%    most functions.
%
%    A new traces object, with no data, can be created as follows:
% 
%          traces = Traces(nTraces,nFrames);
% 
%    This class defines a handle object, which means there is only ever one copy
%    of the data. Changes to any reference to the object will not create a copy;
%    they will change all objects referencing the data. This prevents
%    unnecessary (and potentially dangerous) copying of massive datasets. This
%    behavior also allows subclassing with memory-mapped files.
%
%    To create a copy of an existing object, use the copy() method instead:
%
%          traces2 = copy(traces);
%
%    Some methods used to manipulate the size and shape of the data in place:
%
%       data.subset(idxTraces,idxFrames)   keeps only the specified traces (frames)
%       data.truncate(N)                   truncates all traces to N frames.
%
%    Other methods create a new instance:
%
%       smallData = data.getSubset(idxTraces,idxFrames) create a new instances with only the selected traces (frames)
%       bigData = combine(data,data2,...)   combines all data from multiple instances
% 
%
%    NOTE: for performance reasons, changes to the data fields are not monitored
%    for internal consistency. You can easily create invalid objects without an
%    error. You can even assign and access channels that exist but are not in
%    channelNames (which should mean they are not used), which can confuse other
%    programs. Try to be responsible.
%

% ====== IMPLEMENTATION NOTES =====
% This was not implemented as a Trace class that can be a structure array so
% that the data can be easily memory mapped.
%
% I tried using get/set functions to abstract an underlying cell array of fields
% or a multi-dimensional array, but the overhead is significant and copying is a
% problem in set methods. Overloading subsref/subsagn works and is fast, but
% accessing the traceMetadata fields doesn't act right -- you cannot do this:
%    donor_x = [data.traceMetadata.donor_x]
% This may be solvable with a better-written subsref....
% The advantage of doing that would be a way to iterate over fluorescence
% channels without caring about their names, simplifying some code.
%



% Declaration and initialization of public parameters
properties (Access=public)
    % Data channels (donor,acceptor,fret,etc) must be defined in a subclass
    % implementing this interface. The names of those properties must be the
    % same as those in the channelNames property.
    % donor, acceptor, fret, etc. = [];
    
    time = [];      %time axis in ms or frames
    
    fileMetadata  = struct([]);  %struct array, one per trace.
    traceMetadata = struct([]);

end %end public properties


% Properties calculated as needed, but not actually saved in the class.
properties (Dependent)
    nChannels; %calculated from channelNames/data
    nFrames;   %calculated from time
    sampling;  %time resolution in milliseconds.
end

% Protected properties set in the constructor and never modified.
% nTraces/nFrames could be Dependent and calculated from the size of the data
% channels, but this implementation solves some other problems.
properties( SetAccess=protected, GetAccess=public );
    nTraces = 0;
    channelNames = {};    % short names for each channel (donor, acceptor, etc)
end




methods
    %% ================ CONSTRUCTORS ================ %%
    function this = Traces(varargin)
        
        % Create an empty object.
        if nargin==0, return; end
        
        if nargin==1 && ischar(varargin{1}),
            % Matlab does not allow a constructor to return an object of a
            % different class (not even a subclass).
            error('Cannot load data into an abstract Traces object');
        end
        
        % -------- Copy constructor --------
        if nargin==1 && strcmp(class(varargin{1}),class(this))
            this = copy(varargin{1});
            return;
        elseif nargin==1,
            error('Invalid constructor parameters for class Traces');
        end
        
        % -------- Constructor with defined traces sizes --------
        % NOTE: get methods for nTraces and nFrames check internal consistency.
        % Until the instance is fully built, do not access those methods!!
        assert( nargin>=2 & isscalar(varargin{1}) & isscalar(varargin{2}), ...
                          'Invalid constructor parameters for class Traces' );
        assert( varargin{1}>=0 & varargin{2}>=0, 'Data size must be positive integers' );
        
        this.nTraces = varargin{1};
        nFrames = varargin{2};
        this.time = 1:nFrames;
            
        % Generate IDs for each trace.
        ids = cellfun(@num2str,num2cell(1:varargin{1}),'UniformOutput',false)';
        if numel(ids)>0
            ids = strcat('Trace#',ids);
        end
        
        this.traceMetadata = struct('ids',ids);
        
        % 
        if nargin>=3 && iscell(varargin{3}),
            assert( all(ismember(varargin{3},properties(this))), ...
                    'Invalid channel name' );
            this.channelNames = varargin{3};
        end

        % Allocate space for the data channels.
        for i=1:numel(this.channelNames),
            this.( this.channelNames{i} ) = zeros(varargin{1},nFrames);
        end
    end
    
    
    
    %% ==============  CHANNEL LIST MANIPULATION and CHECKS  ============== %%
    
    % Verify the object is valid and throw errors if not. This can be called in
    % every set method instead of having a distinct test in each one.
    function checkValid(this)
        sz = [this.nTraces this.nFrames];
        valid = 1;
        
        if ~isempty( setdiff(this.channelNames,properties(this)) ),
            error('Invalid channel names');
        end
        
        for c=1:this.nChannels,
            ch = this.channelNames{c};
            valid = valid & all(size(this.(ch))==sz);
        end
        
        valid = valid & numel(this.traceMetadata)==this.nTraces;
        
        if ~valid,
            error('Channel data dimensions not consistent');
        end
    end
    
    % Verify that the RHS of a data channel set operation is valid.
    % By having this in all set methods, we can insure that all data channels
    % have the same size and match the nTraces/nFrames parameters.
    function checkDataInput(this,val,ch)
        assert( size(val,1)==this.nTraces & size(val,2)==this.nFrames, ...
                    'Invalid data size when altering data property' );
        
        if nargin>2,
            assert( ismember(ch,properties(this)), 'Invalid channel name' );

            % It is valid to set to a channel that is not used (and therefore not in
            % channelNames), but it must still be a valid channel.
            if ~ismember(ch,this.channelNames),
                disp(['Adding channel to Traces object:' ch]);
                this.channelNames = [this.channelNames ch];
            end
        end
    end
    
    % Add a data channel to this instance. Must be one of the allowed (but
    % possibly empty) channels. Does nothing if the channel is already there.
    function addChannel(this,chName,val)
        if ~ismember(chName,properties(this))
            error('Invalid channel name');
        end
        
        if ~ismember(chName,this.channelNames),
            this.channelNames = [this.channelNames chName];
        
            if nargin<3,
                val = zeros(this.nTraces,this.nFrames);
            else
                assert( size(val,1)==this.nTraces, size(val,2)==this.nFrames );
            end

            this.(chName) = val;
        end
    end
    
    
    % Make the class act like a struct for convenience with older code.
    function out = isfield(this,fieldname)
        out = ismember( fieldname, properties(this) );
    end
    
    function out = isChannel(this,fieldname)
        out = ismember(fieldname,this.channelNames);
    end
    
    % When altering the list of channel names, insure the list is valid.
    % Ideally we would initialize any variables in here that are new.
    function set.channelNames(this, chNames)
        
        if ~iscell(chNames) || ~isempty(setdiff(chNames,properties(this))), 
            error('Invalid channel names');
        end
        
        this.channelNames = chNames;
    end
    
    
    %% ================  GET/SET METHODS  ================ %%
    
    function M = get.nFrames(this)
        M = numel(this.time);
    end
    
    function set.nFrames(this,N)
        this.truncate(N);
    end
    
    function C = get.nChannels(this)
       C = numel(this.channelNames);
    end
    
    % Get the time resolution in ms (or frames)
    function dt = get.sampling(this)
        if this.nFrames>1,
            dt = diff( this.time(1:2) );
        else
            dt = [];
        end
    end
    
    function out = isempty(this)
        out = this.nTraces==0;
    end
    
    
    
    %% ================  DATA SLICING and COMBINING  ================ %%
    % NOTE: these will be used by inheriting classes, where channels are
    % defined, but they have no meaning here and will do nothing.
    
    
    % Select a subset of traces (and/or frames)
    function this = subset( this, idxTraces, idxFrames )
        assert( nargin>=2, 'Not enough input arguments' );
        
        % Verify input arguments are valid.
        assert( all(idxTraces>=1 & idxTraces<=this.nTraces), 'Invalid trace indexes' );

        if nargin<3,
            idxFrames = 1:this.nFrames;
        else
            assert( all(idxFrames>=1 & idxFrames<=this.nFrames), 'Invalid frame indexes');
        end
        
        % Extract a subset of traces.
        for c=1:this.nChannels,
            ch = this.channelNames{c};
            this.(ch) = this.(ch)(idxTraces,idxFrames);
        end
        
        this.traceMetadata = this.traceMetadata(idxTraces);
        this.time = this.time(idxFrames);
        this.nTraces = numel(idxTraces);
    end %function subset
    
    
    % Truncate all traces up 
    function this = truncate( this, nFrames )
        assert( nargin>1 );
        assert( all(nFrames>=1 & nFrames<=this.nFrames), 'Invalid frame indexes' );
        assert( numel(idxFrames)>1, 'Cannot create empty Traces objects' );

        % Truncate traces in all channels.
        for c=1:this.nChannels,
            this.(this.channelNames{c}) = this.(this.channelNames{c})(:,1:nFrames);
        end
        
        this.time = this.time(1:nFrames);
    end %function truncate
        
    % Extract a subset of traces from the data. This actually creates a copy of
    % the traces object!
    function data = getSubset(this,idx)
        % This creates a shallow copy -- data elements are not duplicated in
        % memory until they are modified during subset.
        data = copy(this);
        data.subset(idx);
    end
    
    
    % Copy all internal data from some other Traces object into this one.
    % This is useful in constructors, where class type must be preserved.
    % Only copies channels that are common to BOTH objects.
    function copyDataFrom(this,source)
        % Copy all data channels in common with other object to this one.
        for c=1:source.nChannels,
            ch = source.channelNames{c};
            
            if this.isfield(ch),
                this.(ch) = source.(ch);
                
                % If this channel is valid, but not used, add it.
                if ~ismember(ch,this.channelNames)
                    this.channelNames = [this.channelNames ch];
                end
            end
        end
        
        % Copy parameters common to all Traces classes.
        this.nTraces = source.nTraces;
        this.traceMetadata = source.traceMetadata;
        this.fileMetadata = source.fileMetadata;
        this.time = source.time;
    end
    
    
    % Combine Traces objects.
    % Adjustments to trace length or incompatible metdata fields are
    % automatically made so the inputs can be combined.
    % Throughout this function 'obj' refers to the current object in the list of
    % all Traces instances to combine. 'i' iterates over these instances.
    function this = combine( varargin )
        assert( nargin>1 );
        
        % Check for options. The user can specify that the end of shorter traces
        % should be extended instead of truncating everything. The padding is
        % zero -- it used to extend the last value.
        mode = 'truncate';
        if nargin>1 && ischar(varargin{end}),
            if ismember(varargin{end},{'extend','truncate'}),
                mode = varargin{end};
            else
                warning('Invalid option for combining traces files');
            end
            
            objects = varargin{1:end-1};
        else
            objects = varargin;
        end
        
        % Save a copy of the current state. Otherwise we will accedentally
        % overwrite the data when new objects are allocated.
        assert(  all( cellfun(@(x) isa(x,'Traces'), varargin) )  );
        this = objects{1};
        
        for i=1:nargin,
            if objects{i}==this,
                objects{i} = copy(this);
            end
        end
        
        % Determine the total number of traces and verify all the Traces
        % objects are compatible
        traceMetadataFields = fieldnames(this.traceMetadata);
        nTracesEach = zeros(nargin,1);
        nFramesEach = zeros(nargin,1);
        
        for i=1:nargin,
            obj = objects{i};
            
            nTracesEach(i) = obj.nTraces;
            nFramesEach(i) = obj.nFrames;
            
            assert( obj.nChannels==this.nChannels & ...
                    all(strcmp(obj.channelNames,this.channelNames)), ...
                    'Traces objects do not have matching channel names' );
                
            % Find the subset of metadata fields in common with all objects.
            traceMetadataFields = intersect( ...
                      traceMetadataFields, fieldnames(obj.traceMetadata) );
        end
        
        % Determine number of frames to save. If the traces are not all the
        % same length, they must be trunctated.
        nFrames = min(nFramesEach);
                
        if min(nFramesEach)~=max(nFramesEach),
            if strcmp(mode,'truncate'),
                disp('Traces are not the same length! Truncating so they match.');
            elseif strcmp(mode,'extend')
                disp('Traces are not the same length! Extending so they match.');
                nFrames = max(nFramesEach);
            end
        end
        this.nTraces = sum(nTracesEach);
        this.time = this.time(1:nFrames);
        
        % Allocate space for new arrays
        for c=1:this.nChannels,
            this.(this.channelNames{c}) = zeros(this.nTraces,nFrames);
        end
                
        % Copy data from each instance into the (now expanded) output.
        traceSoFar = 0;
        
        for i=1:nargin,
            obj = objects{i};
            idx = traceSoFar + (1:obj.nTraces);
            
            % Remove metadatafields that are not in common.
            fieldsToRemove = setdiff( fieldnames(obj.traceMetadata), traceMetadataFields );
    
            if ~isempty(fieldsToRemove),
                obj.traceMetadata = rmfield( obj.traceMetadata, fieldsToRemove );
                disp('Some metadata fields were removed because they were not found in all Traces instances.');
                %text = sprintf( '%s, ',fieldsToRemove{:} );
            end
            
            % Combine metadata
            if i>1,
                this.traceMetadata = [this.traceMetadata obj.traceMetadata];
            end
            
            % Combine data
            for c=1:this.nChannels,
                ch = obj.channelNames{c};
                this.(ch)(idx,1:nFrames) = obj.(ch)(:,1:nFrames);
            end
            traceSoFar = traceSoFar+obj.nTraces;
        end
                
    end %function combine
    
    
    
end %public methods



end %class Traces






