function data = TracesToBNPh5( filename, varargin )
% TracesTohd5  Load fluorescence/FRET data from a .traces format  and 
% converts them to hdf5 file.
% For the attributes, I used the attribute type from BNP-FRET BINNED
% This file makes .traces compatible with BNP-FRET-Binned



% Get input filenames from user if not given.
if nargin<1 || isempty(filename),
    filename = getFiles;
end

if ~iscell(filename),
    filename = {filename};
end

data = struct();

for i=1:numel(filename),
    % Load input traces file
    dataIn = loadTraces( filename{i} );
    
    % Create output filename automatically
    [p,f] = fileparts( filename{i} );

    % Remove the empty spaces int eh file name 
% Original file name 
oldName = f;

% Remove all spaces from the file name
newName = strrep(oldName, ' ', ''); 

    %
    outname = fullfile( p, [newName 'BNP.h5'] );

    % Save channel data to matlab
    data.time = dataIn.time;
    
    for c=1:dataIn.nChannels,
        ch = dataIn.channelNames{c};
        data.(ch) = dataIn.(ch);
    end
    
   % Store intensities and peak positions in output struct
    donor = data.donor(1,:)';
    acceptor = data.acceptor(1,:)';

acceptor_backgr = 0;
donor_backgr = 0;
offs_acceptor= 0;
offs_donor= 0;
var_acceptor = double(var(acceptor)) ;
var_donor = double(var(donor));
%%% acceptor channel information
acceptor_channel = acceptor; % dataset1
acceptor_channel_backgr = acceptor_backgr; % dataset3

%%% donor channel information
donor_channel = donor; % dataset5
donor_channel_backgr = donor_backgr; % dataset6

offset_acceptor = offs_acceptor; % dataset7
offset_donor = offs_donor; % dataset8
variance_acceptor = var_acceptor; % dataset9
variance_donor = var_donor; % dataset10
    
% Create file
if exist(outname, 'file')
    delete(outname);
end

fid = H5F.create(outname, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');

% Acceptor channel (array)
h5create(outname, '/acceptor_channel', size(acceptor_channel,1));
h5write(outname, '/acceptor_channel', acceptor_channel);


% Create scalar dataspace
space_id = H5S.create('H5S_SCALAR');

% Define datatype (example: double)
type_id = H5T.copy('H5T_NATIVE_DOUBLE');

% Create dataset
dset_id = H5D.create(fid, '/acceptor_channel_bg', type_id, space_id, 'H5P_DEFAULT');
% Write scalar value
H5D.write(dset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', acceptor_channel_backgr);

% Close everything
H5D.close(dset_id);
H5S.close(space_id);
H5T.close(type_id);

% Donor channel (array)
h5create(outname, '/donor_channel', size(donor_channel,1));
h5write(outname, '/donor_channel', donor_channel);

% Donor background (SCALAR)
% Create scalar dataspace
space_id = H5S.create('H5S_SCALAR');

% Define datatype (example: double)
type_id = H5T.copy('H5T_NATIVE_DOUBLE');

% Create dataset
dset_id = H5D.create(fid, '/donor_channel_bg', type_id, space_id, 'H5P_DEFAULT');
% Write scalar value
H5D.write(dset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', donor_channel_backgr);

% Close everything
H5D.close(dset_id);
H5S.close(space_id);
H5T.close(type_id);


% Create scalar dataspace
space_id = H5S.create('H5S_SCALAR');

% Define datatype (example: double)
type_id = H5T.copy('H5T_NATIVE_DOUBLE');

% Create dataset
dset_id = H5D.create(fid, '/offset_acceptor', type_id, space_id, 'H5P_DEFAULT');
% Write scalar value
H5D.write(dset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', offset_acceptor);

% Close everything
H5D.close(dset_id);
H5S.close(space_id);
H5T.close(type_id);


% Create scalar dataspace
space_id = H5S.create('H5S_SCALAR');

% Define datatype (example: double)
type_id = H5T.copy('H5T_NATIVE_DOUBLE');

% Create dataset
dset_id = H5D.create(fid, '/offset_donor', type_id, space_id, 'H5P_DEFAULT');
% Write scalar value
H5D.write(dset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', offset_donor);

% Close everything
H5D.close(dset_id);
H5S.close(space_id);
H5T.close(type_id);

% Create scalar dataspace
space_id = H5S.create('H5S_SCALAR');

% Define datatype (example: double)
type_id = H5T.copy('H5T_NATIVE_DOUBLE');

% Create dataset
dset_id = H5D.create(fid, '/variance_acceptor', type_id, space_id, 'H5P_DEFAULT');
% Write scalar value
H5D.write(dset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', variance_acceptor);

% Close everything
H5D.close(dset_id);
H5S.close(space_id);
H5T.close(type_id);

% Create scalar dataspace
space_id = H5S.create('H5S_SCALAR');

% Define datatype (example: double)
type_id = H5T.copy('H5T_NATIVE_DOUBLE');

% Create dataset
dset_id = H5D.create(fid, '/variance_donor', type_id, space_id, 'H5P_DEFAULT');
% Write scalar value
H5D.write(dset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', variance_donor);

% Close everything
H5D.close(dset_id);
H5S.close(space_id);
H5T.close(type_id);
% Close file
H5F.close(fid);
    

 end
 % Create output filename automatically


end %function LoadTracesBinary









