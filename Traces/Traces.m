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
%          traces = Traces(nTraces,nFrames,channelNames);  %optional
% 
%    This class defines a handle object, which means there is only ever one copy
%    of the data. Unlike normal Matlab variables, if you make several copies of
%    a handle object and change one, all copies are updated simultaneously
%    This prevents unnecessary (and potentially dangerous) copying of massive
%    datasets. This behavior might also allow memory-mapped channel data.
%
%    To create a copy of an existing object, use the copy() method instead:
%
%          traces2 = copy(traces);
%
%    Some methods used to manipulate the size and shape of the data in place:
%
%       data.subset(idxTraces)             %keeps only the specified traces
%       data.subset(idxTraces,idxFrames)   %keeps only the specified traces and frames
%       data.truncate(N)                   %truncates all traces to N frames.
%
%    Other methods create a new instance:
%
%       smallData = data.getSubset(idxTraces,idxFrames)  %create a new instances with only the selected traces (frames)
%       bigData   = combine(data,data2,...)              %combines all data from multiple instances
% 
% See also: TracesFret, TracesFret4, TracesAlex, etc.
%

% ====== IMPLEMENTATION NOTES =====
% NOTE: for performance reasons, changes to the data fields are not monitored
% for internal consistency, so it is possible to create Traces objects that are
% invalid if you are careless. Maintaining internal consistency is the
% responsibility of the user. This may change in the future if Matlab's get/set
% function performance improves.
% 
% I came to this implementation after trying many others, including overriding
% structure-array-like syntax, but all other possible ways have significant
% overhead and are not feasible. In this implementation, the channel data are
% transparent property variables that can be manipulated by the use at will.
%
% Overloading subsref/subsagn works and is fast, but accessing the traceMetadata
% fields doesn't act right -- for example you cannot do this:
%    donor_x = [data.traceMetadata.donor_x]
% This may be solvable with a better-written subsref....
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
        
        assert( isscalar(this.nTraces) && isscalar(nFrames), ...
                'Invalid Traces constructor parameters; sizes must be scalar' );
            
        % Generate IDs for each trace.
        ids = strsplit( sprintf('Trace#%d ',1:this.nTraces) );
        this.traceMetadata = struct( 'ids', ids(1:end-1) );
        
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
        
        % If a logical array is given, convert it into a list of indexes
        if islogical(idxTraces),
            idxTraces = find(idxTraces);
        end
        
        if nargin>2 && islogical(idxFrames),
            idxFrames = find(idxFrames);
        end
        
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
    
    
    % Truncate all traces up so that the final length is nFrames
    function this = truncate( this, nFrames )
        assert( nargin>1 && ~isempty(nFrames) && numel(nFrames)==1 && ...
                isnumeric(nFrames) && nFrames<=this.nFrames, 'Invalid final trace length' );
        this.subset(1:this.nTraces,1:nFrames);
    end %function truncate
    
    
    % Extend traces (using the last value) so that the final length is nFrames
    function this = extend( this, nFrames )
        assert( nargin>1 && ~isempty(nFrames) && numel(nFrames)==1 && ...
                isnumeric(nFrames) && nFrames>this.nFrames, 'Invalid final trace length' );

        % Extend all channels with the last datapoint as padding.
        delta = nFrames-this.nFrames;
        
        for c=1:this.nChannels,
            this.(this.channelNames{c}) = [  this.(this.channelNames{c}) ...
                  repmat( this.(this.channelNames{c})(:,end), [1 delta])  ];
        end
        
        % Extend time axis.
        this.time = this.sampling*(1:nFrames);
    end
    
    
    % Extract a subset of traces from the data. This actually creates a copy of
    % the traces object!
    function data = getSubset(this,varargin)
        % This creates a shallow copy -- data elements are not duplicated in
        % memory until they are modified during subset.
        data = copy(this);
        data.subset(varargin{:});
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
                mode = lower( varargin{end} );
            else
                error('Invalid option for combining traces files.');
            end
            
            objects = varargin(1:end-1);
        else
            objects = varargin;
        end
        
        % Create a copy of the first object as a template for output.
        assert(  all( cellfun(@(x) isa(x,'Traces'), objects) )  );
        this = copy(objects{1});
        
        % Determine the total number of traces and verify all the Traces
        % objects are compatible
        traceMetadataFields = fieldnames(this.traceMetadata);
        nTracesEach = zeros(numel(objects),1);
        nFramesEach = zeros(numel(objects),1);
        
        for i=1:numel(objects),
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
        assert( nFrames>0 );
                
        if min(nFramesEach)~=max(nFramesEach),
            if strcmp(mode,'truncate'),
                this.time = this.time(1:nFrames);
            elseif strcmp(mode,'extend')
                nFrames = max(nFramesEach);
                this.time = this.sampling*(1:nFrames);
            end
        end
        this.nTraces = sum(nTracesEach);
        
        % Allocate space for new arrays
        for c=1:this.nChannels,
            this.(this.channelNames{c}) = zeros(this.nTraces,nFrames);
        end
                
        % Copy data from each instance into the (now expanded) output.
        traceSoFar = 0;
        
        for i=1:numel(objects),
            obj = objects{i};
            tm = obj.traceMetadata;
            idx = traceSoFar + (1:obj.nTraces);
            
            % Remove metadatafields that are not in common.
            fieldsToRemove = setdiff( fieldnames(tm), traceMetadataFields );
            tm = rmfield( tm, fieldsToRemove );
            
            % Combine metadata
            if i>1,
                this.traceMetadata = [this.traceMetadata tm];
            end
            
            % Combine data
            for c=1:this.nChannels,
                ch = obj.channelNames{c};
                nf = min(nFrames,obj.nFrames); %non-extended number of frames.
                this.(ch)(idx,1:nf) = obj.(ch)(:,1:nf);
                
                % If necessary, pad with last data value in every trace.
                if nFrames>obj.nFrames
                    this.(ch)(idx,(nf+1):nFrames) = repmat( obj.(ch)(:,nf), [1 nFrames-obj.nFrames] );
                end
            end
            traceSoFar = traceSoFar+obj.nTraces;
        end
                
    end %function combine
    
    
    
end %public methods



end %class Traces






