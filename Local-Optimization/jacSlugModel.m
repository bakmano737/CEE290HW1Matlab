function J = jacSlugModel(pars,t,Q,d)
% Calculate the Jacobian of the slugmodel at current parameter values

% Define the storage and transmissivity
S = pars(1);
T = pars(2);

% Now derive Jacobian at each "t" value

% Derivative with respect to "S"
%J(:,1) = ( ( (-1/16 *Q / pi / T^2) ./ t.^2 * d^2 ) .* exp( - 1/4 * d^2 * S / T ./t) )';
J(:,1) =  ( ( ( 1/16 *Q / pi / T^2) ./ t.^2 * d^2 ) .* exp( - 1/4 * d^2 * S / T ./t) )';

% Derivative with respect to "T"
%J(:,2) = ( ( -1/4 * Q / pi / T^2 ./t ) .* exp ( -1/4 * d^2 * S / T ./t) + (1/16 * Q / pi / T^3 ./t.^2 * d^2 * S).* exp( -1/4 * d^2 * S / T ./t) )';
J(:,2) =  ( (  1/4 * Q / pi / T^2 ./t ) .* exp ( -1/4 * d^2 * S / T ./t) - (1/16 * Q / pi / T^3 ./t.^2 * d^2 * S).* exp( -1/4 * d^2 * S / T ./t) )';