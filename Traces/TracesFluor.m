classdef TracesFluor < Traces
% TracesFret: multi-channel fluorescence traces
% 
%    Traces class with generic fluorescence channels and no derived data
%    channels (such as FRET efficiency). Data channels have simple names of
%    ch1, ch2, ch3, etc. For now, this is limited to 4 channels at most.
%
%    These data channels may be empty if any are not used. When querying the
%    trace data, you can iterate over all channels using 'channelNames' or over
%    a specific type of data using idxFret and idxFluor.
%
%    A new traces object, with all data filled with zeros, can be created:
% 
%          traces = Traces(nTraces,nFrames,nChannels);
%
%    See Traces.m for more information on Traces objects.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Declaration and initialization of public parameters
properties (Access=public)
    ch1 = [];
    ch2 = [];
    ch3 = [];
    ch4 = [];
end %end public properties




methods
    
    
    %% ================ CONSTRUCTORS ================ %%
    function this = TracesFluor(varargin)
        % NOTE: I did not try to implement a copy constructor. There may be
        % limitations in the language that prevents this.
        assert(nargin>=2);
        nTraces = varargin{1};
        nFrames = varargin{2};
        
        if nargin<3
            nCh = 4;
        else
            nCh = varargin{3};
        end
        assert( nCh>=1 | nCh<=4, 'Invalid number of fluorescence channels. Must be 1..4' );
        
        % Generate generic channel names: ch1, ch2, etc.
        chNames = cell(nCh,1);
        for c=1:nCh,
            chNames{c} = sptrintf('ch%d',c);
        end
        
        % Call superclass constructor to set up the object
        this = this@Traces(nTraces,nFrames,chNames);
    end
    
    
     %% ================  GET/SET METHODS  ================ %%
     
    % Calculate total fluorescence intensity
    function T = total(this)
        T = zeros(this.nTraces,this.nFrames);
        for c=1:this.nChannels,
            T = T + this.( this.channelNames{c} );
        end
    end
     
    
end %public methods

end %class Traces






