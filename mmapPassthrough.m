classdef mmapPassthrough
% mmapPassthrough  class to transparently memory map data.
% This class simple acts as a wrapper around a memory mapped variable.
% Any references to elements or properties of the class are passed directly to
% the memory mapped variable (Data element in the memmapfile object).
% To the user, this will seem as if the object is the memory mapped variable
% itself.
%
% Example usage:
%   mmap = memmapfile( filename, 'Format', {'uint16',[3,5],'x'} );
%   passthrough = mmapPassthrough( mmap, options );
%
% Several options can be specified:
%  .forwardToVar: if there is only one variable in the 'Format' of the memory
%       map, it can be treated as the data to be wrapped. Give the name of
%       the variable ('x').
%  .rowMajor: set as true if the underlying data in the file was written in
%       row-major order. When creating the memory map, you MUST reverse the
%       order of width and height for the map size to account for this.
%
    
% Private properties.
properties  (SetAccess='private', GetAccess='private')
    mmap = [];
    rowMajor = false;  %transpose subscripts for row-major file data.
    forwardToVar = []; %forward references to this member variable of Data.
end %properties


% Public system methods 
methods

    %---------       Constructor        ---------%
    function obj = mmapPassthrough( MMAP, options )
        obj.mmap = MMAP;
        
        % Process optional arguments
        if nargin>1 && isfield(options,'rowMajor'),
            obj.rowMajor = options.rowMajor;
        end
        if nargin>1 && isfield(options,'forwardToVar'),
            obj.forwardToVar = options.forwardToVar;
        end
        
    end %constructor
    
    
    
    %---------       Direct data element reference        ---------%
    
    % Any reference to properties of this class is passed straight through
    % to the underlying object (data).
    function B = subsref(obj,S)
        
        % If the file data are in row-major order, re-order the subscripts of
        % the first two dimensions to translate from column-major (MATLAB).
        if obj.rowMajor && S(1).type(1)=='(',
            S(1).subs([1,2]) = S(1).subs([2,1]);
        end
        
        %
        if isempty(obj.forwardToVar),
            B = subsref( obj.mmap.Data, S ); 
        else
            B = subsref( obj.mmap.Data.(obj.forwardToVar), S );
        end
        
        % If accessing in row-major order, we have no choice but to transpose
        % the data for access in MATLAB. A copy is basically forced.
        if obj.rowMajor,
            B = ndTranspose(B);
        end
    end
    
    function obj = subsasgn(obj,S,val)        
        % If the file data are in row-major order, re-order the subscripts of
        % the first two dimensions to translate from column-major (MATLAB).
        if obj.rowMajor && S(1).type(1)=='(',
            S(1).subs([1,2]) = S(1).subs([2,1]);
        end
        
        % Build the subscript series to get to the mmap'd data.
        S(2:numel(S)+1) = S(:);
        S(1).type = '.';
        S(1).subs = 'Data';
        
        % 
        if ~isempty(obj.forwardToVar),
            S(2:numel(S)+1) = S(:);
            S(1).type = '.';
            S(1).subs = obj.forwardToVar;
        end
        
        % If accessing in row-major order, a transpose is necessary.
        % WARNING: Untested!Direct data element reference
        if obj.rowMajor,
            val = ndTranspose(val);
        end
        
        % Pass the assignment on to memmapfile.
        subsasgn( obj.mmap, S, val );
    end
    
    
    
    %---------       Data size query methods        ---------%

    function varargout = size( obj, dim )
        
        if isempty(obj.forwardToVar),
            S = size( obj.mmap.Data ); 
        else
            S = size( obj.mmap.Data.(obj.forwardToVar) );
        end
        
        
        if nargin>1,
            assert( nargout<=1, 'Too many output arguments' );
            varargout{1} = S(dim);
            return;
        end
        
        if nargout<=1,
            varargout{1} = S;
        elseif nargout<=numel(S),
            varargout = num2cell( S );
        else
            error( 'Too many output arguments' );
        end
    end
    
    function N = numel( obj, varargin )
        error( nargchk(1,1,nargin) );
        
        N = prod( size(obj) );
    end
    
    function N = end(obj,k,n)
        assert( n==numel(size(obj)) );
        
        N = size(obj,k);
    end
    
end %public methods



end %classdef



function Bnew = ndTranspose( B )
% Transposes the first two dimensions (width/height) of the matrix B, while
% leaving the order of the other elements the same.
% WARNING: does not support more than 3 dimensions!

    Bnew = zeros( size(B,2), size(B,1), size(B,3), class(B) );

    for i=1:size(B,3),
        Bnew(:,:,i) = B(:,:,i)';
    end
end
