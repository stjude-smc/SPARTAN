classdef TracesFret4 < TracesFret
% TracesFret: multi-channel FRET traces
% 
%    Traces class with extra properties for additional data channels, including
%    donor2, acceptor2, and fret2 for three-color FRET and two-pair FRET
%    experiments, and factor/factor2 channels for miscellaneous fluorescence
%    channels (for example for dye-labeled factor binding to the ribosome).
%
%    These data channels may be empty if any are not used. When querying the
%    trace data, you can iterate over all channels using 'channelNames' or over
%    a specific type of data using idxFret and idxFluor.
%
%    A new traces object, with all data filled with zeros, can be created as
%    follows, where the default channels are donor, acceptor, and fret.
% 
%          traces = Traces(nTraces,nFrames);
%
%    If only a subset of fields is valid, add an argument for channel names:
%
%          traces = Traces(nTraces,nFrames,channelNames);
%
%
%    See Traces.m and TracesFret.m for more information on Traces objects.
% 



% Declaration and initialization of public parameters
properties (Access=public)
    % Other fields for multi-color fret. These may be empty if they are not
    % used. To know what fields are "alive", check the channelNames list or find
    % which of these fields are not empty.
    donor2    = [];
    acceptor2 = [];
    fret2     = [];
    factor    = [];
    factor2   = [];
end %end public properties


properties (Dependent)
    idxFactor;  %indexes (in channelNames) of factor-binding channels.
end



methods
    %% ================ CONSTRUCTORS ================ %%
    function this = TracesFret4(varargin)
        % Constructors cannot create an object of a different class than their
        % own. Instead the data are copied into an instance of this class.
        if nargin>=1 && (ischar(varargin{1}) || isa(varargin{1},'Traces')),
            args = {0,0,{}};
        else
            args = varargin;
        end
        
        % Call superclass constructor to set up the object.
        this = this@TracesFret( args{:} );
                
        % -------- Load traces from file --------
        if nargin==1 && ischar(varargin{1}),
            this.copyDataFrom( loadTraces(varargin{1}) );
            
        % -------- Copy Constructors --------
        % This does not allow initializing TracesFret4 object with TracesFret
        % data, which seems like a valid thing to do! Doing so may require
        % copying property fields.
        elseif strcmp(class(varargin{1}),class(this))
            this = copy(varargin{1});
            
        elseif nargin==1 && isa(varargin{1},'Traces')
            this.copyDataFrom(varargin{1});
        end
    end
    
    
     %% ================  GET/SET METHODS  ================ %%
     
     function idx = get.idxFactor(this)
        idx = find(  cellfun( @(x) ~isempty(strfind(x,'factor')), this.channelNames )  );
     end
     
    
    
end %public methods

end %class Traces






