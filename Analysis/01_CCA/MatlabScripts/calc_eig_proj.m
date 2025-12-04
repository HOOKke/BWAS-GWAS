function [New_Eigvec] = calc_eig_proj(NewData,OrigData,Eigvec,Eigval,varargin)

assert(isequal(size(NewData),size(OrigData)),...
    'New and old data must have same dimension.')
p = inputParser;
addParameter(p,'thresh',1e-8)
parse(p,varargin{:})

thresh = p.Results.thresh;

if isvector(Eigval)
    Eigval = diag(Eigval);
end

Xnew = NewData;
X = OrigData;
U = Eigvec;
D = Eigval;     % note: should be eigendecomp of XX^*, NOT of X

D = D.*(D >= thresh);
keep = sum(D) > 0;
D = D(:,keep);

Xp_old = pinv(X);
Xp_new = pinv(Xnew);

% S = sqrt(D);
% V_orig = Xp_old*U*S;

% New_Eigvec = Xp_new'*V_orig*S;
New_Eigvec = Xp_new'*Xp_old*U*D;
end