function Q   = sampler_full_rate(P,escrate,params)

if ~isrow(escrate)
    keyboard
end

% keyboard
Q                    = -diag(escrate) + P.*escrate';
