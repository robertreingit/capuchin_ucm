function ucm_res = generic_UCM(Jacobian, params, take_sqrt)
% generic Uncontrolled manifold function
% calculates the uncontrolled manifold as described in
%   Scholz, JP & Schöner, G (1999). The uncontrolled manifold concept:
%       Indentifying control variables for a functional task,
%       Exp Brain Res, 126: 289-306.
% In the orginal paper, the square-root of the squared and summed variance 
% was taken before further analysis, however, in more recent analysis this.
% is not the case any more. Therefore, the last parameter specifies
% whether the square root of the raw values should be taken.
% 
% INPUT:
%   Jacobian: Jacobian matrix
%   params: input variables
%   take_sqrt: Boolean, take the square root of variance {default: true}
% OUTPUT:
%   ucm_res: result structure
%       .ucm_par:
%       .ucm_ort:
%       .var_par_ind:
%       .var_ort_ind:
% SIDEEFFECTS:
%  None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
    take_sqrt = true;
end

[NO_TRIALS, NO_PARAMETERS] = size(params);
d = size(Jacobian, 1);

% calculate null space
ucm = null(Jacobian);

ucm_par_total = 0;
ucm_ort_total = 0;
ucm_total_total = 0;

Omega_p = zeros(1, NO_PARAMETERS);
Omega_o = zeros(1, NO_PARAMETERS);

Omega_bar = mean(params);
for trial = 1:NO_TRIALS
    
    dOmega = params(trial,:) - Omega_bar;
    % vector parallel to UCM.
    % Must be column vector for latter stuff to be commensurate.
    Omega_p_i = pV2Sub(dOmega, ucm)';
    % vector orthogonal to UCM.
    Omega_o_i = dOmega - Omega_p_i;
    
    % calculate variation per parameter
    Omega_p = Omega_p + Omega_p_i.^2;
    Omega_o = Omega_o + Omega_o_i.^2;
    
    % accumulate squared deviations.
    ucm_par_total = ucm_par_total + sum(Omega_p_i.^2);
    ucm_ort_total = ucm_ort_total + sum(Omega_o_i.^2);
    ucm_total_total = ucm_total_total + sum(dOmega.^2);
    
end


% check whether sqrt should be applied
if take_sqrt
    collect = @(x) sqrt(x);
else % just identity function
    collect = @(x) x;
end

% calculate result parameters
ucm_par = collect(ucm_par_total / ((NO_PARAMETERS - d) * NO_TRIALS));
ucm_ort = collect(ucm_ort_total / (d * NO_TRIALS));
ucm_tot = collect(ucm_total_total / (NO_PARAMETERS * NO_TRIALS));
Omega_p = collect(Omega_p / NO_TRIALS);
Omega_o = collect(Omega_o / NO_TRIALS);

% build result structure
ucm_res.par_total = ucm_par;
ucm_res.ort_total = ucm_ort;
ucm_res.tot_total = ucm_tot;
ucm_res.var_par_ind = Omega_p;
ucm_res.var_ort_ind = Omega_o;
