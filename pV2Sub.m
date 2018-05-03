function pvec = pV2Sub(v, ONB)
% pV2Sub(v,ONB)
% Projects the vector v onto the subspace ONB
% INPUT:
%   v:      vector to be projected
%   ONB:    ORTHONORMAL(!!!) basis vectors, routine doesn't check
% OUTPUT:
% pvec = projected vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isrow(v)
    v = v';
end

pvec = ONB*(ONB'*v);
