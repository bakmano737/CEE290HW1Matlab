function h_pred = slugmodel(pars,t,Q,d)
% Run the forward model
%

% Define the storage and transmissivity
S = pars(1);
T = pars(2);

% Predict "h" using the slug model equation (see reader)
h_pred = ( Q./ (4 * pi * T * t) ) .* exp ( ( -d^2 * S ) ./ ( 4 * T * t) );