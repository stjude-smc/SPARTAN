classdef mmapPassthrough
% mmapPassthrough:
% This class simple acts as a wrapper around a memory mapped variable.
% Any references to elements or properties of the class are passed directly to
% the memory mapped variable (Data element in the memmapfile object).
% To the user, this will seem as if the object is the memory mapped variable
% itself.
    
% Private properties.
properties  (SetAccess='private', GetAccess='private')
    mmap = [];
end %properties


% Public system methods 
methods

    % Constructor
    function obj = mmapPassthrough( MMAP )
        obj.mmap = MMAP;
    end %constructor
    
    % Any reference to properties of this class is passed straight through
    % to the underlying object (data).
    function B = subsref(obj,S)
        B = subsref( obj.mmap.Data, S );
    end
    
    function obj = subsasgn(obj,S,val)
        % Build the subscript series to get to the mmap'd data.
        Sc(2:numel(S)+1) = S(:);
        Sc(1).type = '.';
        Sc(1).subs = 'Data';
        disp(Sc);
        
        % Pass the assignment on to memmapfile.
        subsasgn( obj.mmap, Sc, val );
    end
    
end %public methods



end %classdef


