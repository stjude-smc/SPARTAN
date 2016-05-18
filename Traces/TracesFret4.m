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

%   Copyright 2007-2015 Cornell University All Rights Reserved.


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


properties (Constant,Hidden=true)
    % Defines the possible methods for FRET calculcation. The descriptions are
    % used when asking the user to specify what calculation to be used and what
    % assumptions can be made.
    % Hidden just to keep down the clutter when displaying objects.
    fretGeometryNames = {'tandem3','independent3','acceptor/total'};
    fretGeometryDesc  = {'Tandem 3-color (D->A1->A2)','Independent 3-color (D->A1, D->A2)','Acceptor/Total'};
end

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
        elseif nargin==1 && strcmp(class(varargin{1}),class(this))
            this = copy(varargin{1});
            
        elseif nargin==1 && isa(varargin{1},'Traces')
            this.copyDataFrom(varargin{1});
        end
    end
    
    
     %% ================  GET/SET METHODS  ================ %%
     
     function idx = get.idxFactor(this)
        idx = find(  cellfun( @(x) ~isempty(strfind(x,'factor')), this.channelNames )  );
     end
     
    
     function changed = verifyFretGeometry(this)
     % Ensure that the method for calculating FRET is specified for 3/4-color
     % FRET. If not or if invalid, ask the user to specify.
        changed = false;
     
        if ~isChannel(this,'acceptor2'), return; end
     
        if ~isfield(this.fileMetadata,'fretGeometry') || ...
           ~ismember(this.fileMetadata.fretGeometry,this.fretGeometryNames),

            fretGeometry = 'acceptor/total';  %default

            % Handle legacy tag that specifies tandem3 calculation.
            if isfield(this.fileMetadata,'isTandem3') && this.fileMetadata.isTandem3,
                fretGeometry = 'tandem3';

            % Otherwise, we have to ask the user what to do.
            else
                [sel,ok] = listdlg('PromptString','What is the 3-color FRET geometry in this experiment?', ...
                               'SelectionMode','single','ListSize',[300 120], ...
                               'ListString',TracesFret4.fretGeometryDesc, ...
                               'InitialValue',3 );
                if ok,
                    fretGeometry = TracesFret4.fretGeometryNames{sel};
                end
            end

            this.fileMetadata(1).fretGeometry = fretGeometry;
            changed = true;
        end
     end
     
     
     function name = fretAxisLabel(this)
     % Y-axis label of FRET data for display.
     % FIXME: does not handle single-pair with factor channel correctly.
        if ~isfield(this.fileMetadata,'fretGeometry') || numel(this.idxFret)<2,
            name = fretAxisLabel@TracesFret(this);
        elseif strcmpi(this.fileMetadata.fretGeometry,'acceptor/total')
            name = 'Acceptor/Total';
        else
            name = [this.fileMetadata.fretGeometry(1:end-1) ' FRET'];
            name(1) = upper(name(1));
        end
     end
     
     
    %% ================ DATA MANIPULATION ================ %%
    
    function this = recalculateFret( this, varargin )
    % Recalculate FRET efficiencies (all fields) using fluorescence.
    %
    %   data.recalculateFret() recalculates FRET for all traces.
    %   data.recalculateFret(IDX) only alters the subset of traces IDX.
    %   data.recalculateFret(IDX,THRESH) uses the supplied thresholds.
    % 
    %   fileMetadata.fretGeometry specifies the FRET calculation method.
    %   If not available, the user is prompted.
    %     'tandem3':        D->A1, A1->A2. Assumes no D->A2 FRET.
    %     'independent3':   D->A1, D->A2.  Assumes no A1->A2 FRET.
    %     'independent4':   D->A1, D2->A2. Assumes independent FRET pairs.
    %     'acceptor/total': Fraction of total intensity from each acceptor.
    %
    %   fileMetadata.zeroMethod specifies the method to use for setting FRET to
    %   zero when the donor is dark: 'skm' or 'threshold'.
    
        %------------------  Actually calculate FRET  ------------------%
        
        % If there is only one FRET channel, just use the subclass method.
        % This can happen if there are factor channels, but only one FRET pair.
        if ~isChannel(this,'acceptor2'),
            recalculateFret@TracesFret( this, varargin{:} );
            return;
        end
        assert( this.isChannel('fret') && this.isChannel('fret2') );
        
        % Ask the user for the FRET calculation method if not in metadata.
        this.verifyFretGeometry();
        
        % Select the subset of traces to be corrected.
        if nargin>=2,
            idx = varargin{1};
            if isempty(idx), return; end
            assert( isvector(idx) );
            selData = getSubset(this,idx);
        else
            selData = this;
            idx = 1:this.nTraces;
        end
        idx = to_row(idx);
        
        % Calculate FRET for each trace.
        total = selData.total;
        
        switch this.fileMetadata.fretGeometry,
        case 'tandem3',
            % D1->A1->A2 FRET. Assumes no D->A2 FRET.
            acc = selData.acceptor + selData.acceptor2; %total acceptor intensity
            newfret  = acc./total;
            newfret2 = selData.acceptor2./acc;
        case 'independent3',
            % D1->A1, D1->A2 FRET. Assumes no A1->A2 FRET.
            newfret  = selData.acceptor ./(selData.acceptor +selData.donor);
            newfret2 = selData.acceptor2./(selData.acceptor2+selData.donor);
        case 'independent4',
            % D1->A1, D2->A2 FRET. Assumes the two pairs are fully independent.
            %newfret  = selData.acceptor ./(selData.acceptor +selData.donor);
            %newfret2 = selData.acceptor2./(selData.acceptor2+selData.donor2);
            error('Four-color FRET not supported.');
        case 'acceptor/total',
            % Fraction of total intensity in each channel (Not FRET).
            newfret  = selData.acceptor./total;
            newfret2 = selData.acceptor2./total;
        end
        
        % Remove any NaN values. Sometimes happens with low SNR traces.
        newfret( isnan(newfret) ) = 0;
        newfret2( isnan(newfret2) ) = 0;
        
        
        %-------------  Set FRET to zero when the donor is dark  ------------%
        
        % Set FRET to zero when the donor is dark (total intensity at baseline).
        if ~isfield(this.fileMetadata,'zeroMethod')
            this.fileMetadata.zeroMethod = 'threshold';
        end
        
        switch lower(this.fileMetadata.zeroMethod)
            case 'off'
                alive = thresholdTotal( total, zeros(this.nTraces,1) );
            case 'threshold'
                alive = thresholdTotal( total, varargin{2:end} );
            case 'skm'
                alive = skmTotal( total, varargin{2:end} );
            otherwise
                warning('Unknown value for fileMetadata.zeroMethod. Defaulting to threshold method.');
                this.fileMetadata.zeroMethod = 'threshold';
                alive = thresholdTotal( total );
        end
        
        newfret(~alive) = 0;
        newfret2(~alive) = 0;
        
        % For tandem3 configuration, set A1->A2 to zero if the D->A1 efficiency
        % is too low for an accurate value.
        if strcmp(this.fileMetadata.fretGeometry,'tandem3'),
            constants = cascadeConstants;
            
            % Blinking
            newfret2( newfret<0.2 ) = 0;
            
            % Photobleaching
            for i=1:size(newfret,1),
                fretRange = rleFilter( newfret(i,:)>=0.25, constants.rle_min );
                fret2_end = find(fretRange,1,'last');
                if ~isempty(fret2_end)
                    newfret2(i, fret2_end:end ) = 0;
                end
            end
        end
        
        % Save the recalculated subset of traces.
        this.fret(idx,:) = newfret;
        this.fret2(idx,:) = newfret2;
    end
end %public methods

end %class Traces






