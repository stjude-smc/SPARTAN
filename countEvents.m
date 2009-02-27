function [ntrans,tot_time] = countEvents(dwtfilename, edges, nsegments)
%COUNTEVENTS  Gets number of transitions of each class per trace
% 
%   [NT,TIME] = COUNTEVENTS(DWT, EDGES)
%   Uses the idealization in DWT to create a NxM matrix (in NT) of the
%   number of transitions of each class, where N is the number of traces
%   and M is the number of possible transitions (EDGES).  Also returns the
%   total time (in sec) in states above 1 (blinking).
%
%   If EDGES is not specified, NT is simply the total number of
%   transitions.
% 

% TODO: ....
%
% NOTE: dwt lists states as 0-based


if ~exist('nsegments','var'), nsegments=1; end



% Load idealization data file
fid = fopen(dwtfilename);


tot_time = 0;
% ntrans = 0;


% No model givenm just sum up total number of transitions
if ~exist('edges','var')
    ntrans = 0;
    
    while 1,
        textscan(fid, 'Segment: %d Dwells: %d Sampling(ms): %d Start(ms): %d %*[^\n]');

        data = textscan(fid, '%d%d');
        states = data{1};
        times  = data{2};

        if numel(times) == 0,  break;  end  %end of file
        
        ntrans = ntrans + numel(times);
        tot_time = tot_time + sum( times(states>0) );
    end
    
% Model given - produce transition matrix
else
    edges  = edges-1;
    nedges = length(edges);
    
    ntrans = zeros( nsegments, nedges );

    i = 1;
    while 1,
        % Load state data
        textscan(fid, 'Segment: %d Dwells: %d Sampling(ms): %d Start(ms): %d %*[^\n]');

        data = textscan(fid, '%d%d');
        states = data{1};  %States are 0 based, but so are the edges we define
        times  = data{2};

        if numel(states) == 0,  break;  end  %end of file

        tot_time = tot_time + sum( times(states>0) );

        % Count transitions on each edges
        for j=1:nedges,
            ntrans(i,j) = sum( ...
                states(1:end-1)==edges(j,1) & states(2:end)==edges(j,2)  );
        end

        i = i+1;
    end

end %if edges given


fclose(fid);

tot_time = tot_time/1000; %in seconds

end  % function CountEvents


