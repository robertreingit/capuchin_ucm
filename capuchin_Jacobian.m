function J = capuchin_Jacobian(segment_lengths, parameters)
% capuchin_Jacobian
% Calculates the Jacobian for the capuchin-8-segment model.
% INPUT:
%   segment_lengths:    Length of body segments
%   parameters:     1 x 10 vector of parameters
% OUTPUT:
%   J: Jacobian matrix e R^2x10
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

angles = parameters;
if iscolumn(angles)
    angles = angles';
end

cangles = cumsum(angles);
 tmp_mat = [
    -sin(cangles) .* segment_lengths; ...
    cos(cangles) .* segment_lengths ...
    ];
J_tmp = tmp_mat * tril(ones(8));
J = [J_tmp(1,:); ...
    J_tmp(2,:)];
