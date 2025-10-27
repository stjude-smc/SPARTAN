function spartanTelemetry()
% spartanTelemetry  Ping server to collect basic usage statistics.
%

%   Copyright 2025 All Rights Reserved. 


api_key = 'phc_cpOjVvNuHOKGji8EFCs9hGCKToPLYCe960qFAamlI20';
endpoint = 'https://us.i.posthog.com/i/v0/e/';


% Send telemetry at most once per day.
persistent checktime;
if ~isempty(checktime) && datetime('now')<checktime
    return;
end
checktime = datetime('now') + days(1);


% Check for permanent opt out of telementry.
if ~ispref('SPARTAN','out_out')
    setpref('SPARTAN', 'out_out', false);
end
if getpref('SPARTAN', 'out_out')==true
    return;
end


% Get or generate unique installation ID
if ispref('SPARTAN', 'distinct_id')
    distinctID = getpref('SPARTAN', 'distinct_id');
else
    distinctID = char(java.util.UUID.randomUUID);
    setpref('SPARTAN', 'distinct_id', distinctID);
end

% Send payload to PostHog, failing silently
payload = struct( ...
    'api_key', api_key, ...
    'event', 'SPARTAN_started', ...
    'distinct_id', distinctID, ...
    'properties', struct( ...
        'app_version', cascadeConstants('version'), ...
        'platform', computer, ...
        'matlab_version', version, ...
        'isdeployed', isdeployed ...
    ) ...
);

options = weboptions('MediaType','application/json', 'Timeout', 2);
try
    response = webwrite(endpoint, payload, options);
    disp(response);
catch ME
    fprintf('Telemetry send failed: %s\n', ME.message);
end



end  %function




