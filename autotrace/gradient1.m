function y = gradient1(x)
% Gradient of a vector signal (optimized version of gradient)

% note we assume the vector is more than 2 frames long?

y = [  x(2)-x(1)  0.5*(x(3:end)-x(1:end-2))  x(end)-x(end-1)  ];
% y = [  x(:,2)-x(:,1)  0.5*(x(:,3:end)-x(:,1:end-2))  x(:,end)-x(:,end-1)  ];

end