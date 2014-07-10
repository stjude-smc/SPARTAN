function [picks,rejects] = splitPeaks2( total_picks, fieldSize, quadrant, tform )
% Predict peak locations for a specific fluorescence field given the
% locations of the peaks in the total fluorescence intensity image.
% This involves translating the image into the target quadrant and 
% applying a transformation for software alignment.
%
%     picks  = splitPeaks( total_picks, fieldSize, quadrant, [tform] )
%
%     total_picks = peak locations in total intensity imge ([x1,y1;x2,y2;...])
%
%     fieldSize = [nrow,ncol] size of field-of-view.
%
%     quadrant = defines which area of the full field-of-view to
%         translate the peak locations to. For 2-color, 1=left, 2=right.
%         For four-color, 1=UL,2=UR,3=LL,4=LR.
%
%     tform  = transformation (from maketform) to apply to the peak
%         locations to correct for misalignment. This may place some peaks
%         outside of the range! Such peaks are given in "rejects" output.
%         This parameter is optional.
%
%     picks   = peak locations for all fields, with the locations for each
%         each field interleaved (ch1_mol1,ch2_mol1,ch3_mol1,ch1_mol2,ch2_mol2,...).
%         First column is x, second column is y coordinate.
% 
%     rejects = list of peak locations (by index into the total_picks
%         input array) that are projected outside the range of the field
%         because of the software alignment applied. These
%

assert( nargin>=3 & size(total_picks,2)==2 );


% Define the size of the field-of-view.
nrow = fieldSize(1);
if numel(fieldSize)>1,
    ncol = fieldSize(2);
else
    ncol = fieldSize(1);
end


if nargin>3 && ~isempty(tform),
    % Apply transformation to the peak locations.
    % We have to first center the peak locations at (0,0) so the rotation is
    % about the center of the image, then put it back afterward.
    % Peak (maxima) locations must be integers, so they are rounded.
%     picks = tformfwd(  tform,  [total_picks(:,1)-(ncol/2) total_picks(:,2)-(nrow/2)]  );
%     picks = round( [picks(:,1)+(ncol/2) picks(:,2)+(nrow/2) ] );
    
    % Centering isn't required (it actually breaks it) for a real tform.
    picks = round( tformfwd( tform, total_picks ) );
    
    % Mark any peaks that now fall outside the field limits.
    rejects = picks(:,1)<3      | picks(:,2)<3       | ...
              picks(:,1)>ncol-2 | picks(:,2)>nrow-2;
else
    picks = total_picks;
    rejects = false( size(picks,1),1 );
end


% Translate into target field
switch quadrant
    case 1,
        % UL. nothing to do.
        
    case 2,
        picks(:,1) = picks(:,1) + ncol; % UR x

    case 3,
        picks(:,2) = picks(:,2) + nrow; % LL y
        
    case 4,
        picks(:,1) = picks(:,1) + ncol; % LR x
        picks(:,2) = picks(:,2) + nrow; % LR y
        
    otherwise
        error('Invalid quadrant. Must be 1-4.');
end




end  %function splitPeaks
