function y = std1(x)
% standard deviation of a vector signal (optimized version of std)

n = numel(x);

if n>1
    denom = n-1;
else
    denom = n;
end

x = x-(sum(x)./n);
y = sum(x.^2) ./ denom;
y = sqrt(y);

% This is faster, but apparently gives slightly different results.
% inlining sqrt() is somehow different (?!)
% y = sqrt(  sum((x-(sum(x)./n)).^2) ./ denom  );

end


