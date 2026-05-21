function gettraces_fromTwoStacks()
%GETTRACES_FROMTWOSTACKS  Pick two channel TIFFs, merge L-R, open gettraces.
%
%   Requires combineTwoChannelStacks_tile.m on the path.
%   For auto-open of the merged movie, complete Phase 3 in gettraces_gui.m.

filter = {'*.tif;*.tiff', 'TIFF stacks (*.tif, *.tiff)'};

% --- Channel 1 (left) ---
[f1, p1] = uigetfile(filter, 'Two-channel merge: select channel 1 (left)');
if isequal(f1, 0), return; end
path1 = fullfile(p1, f1);

% --- Channel 2 (right) ---
[f2, p2] = uigetfile(filter, 'Two-channel merge: select channel 2 (right)');
if isequal(f2, 0), return; end
path2 = fullfile(p2, f2);

% --- Output path ---
[~, base1, ext1] = fileparts(f1);
defaultName = [base1 '_mergedLR' ext1];
[fOut, pOut] = uiputfile({'*.tif;*.tiff', 'TIFF (*.tif, *.tiff)'}, ...
    'Save merged side-by-side movie as', fullfile(p1, defaultName));
if isequal(fOut, 0), return; end
outPath = fullfile(pOut, fOut);

% --- Merge ---
h = waitbar(0, 'Merging TIFF stacks...');
cleanup = onCleanup(@() closeWaitbar(h)); %#ok<NASGU>
try
    combineTwoChannelStacks_tile(path1, path2, outPath);
catch ME
    closeWaitbar(h);
    errordlg(ME.message, 'Merge failed');
    rethrow(ME);
end
closeWaitbar(h);

% --- Open gettraces ---
if exist('gettraces_gui', 'file') == 2
  %pass merged path so OpeningFcn calls OpenStk
  gettraces_gui(outPath);
else
  gettraces;
  warndlg(['Merged file written to:' newline outPath newline newline ...
      'Use Open Movie and select this file.'], 'Merge complete');
end

end

function closeWaitbar(h)
if ~isempty(h) && ishandle(h)
    close(h);
end
end