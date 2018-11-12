function [X_opt, iter] = capped_tnn_admm(S, X, Y, M, omega, opts ,theta)
%————————theta——————————%
% min_X ||S||_* / theta - < U , X >  s.t.  X - S = 0    
% Input:
%       S  X  Y  M  omega  opts  theta

% Output:
%       X_opt   -    recovered tensor
%       iter    -    numebr of iterations
%--------------------------------------------------------------------------
DISPLAY_EVERY = 10;
rho = 1.05; 
mu = 1e-3; 
max_mu = 1e10; 
admm_tol = 1e-4;
admm_iter = 200;%500？？？？

if ~exist('opts', 'var'),   opts = [];  end
if isfield(opts, 'rho');        rho = opts.rho;              end
if isfield(opts, 'mu');         mu = opts.mu;                end
if isfield(opts, 'max_mu');     max_mu = opts.max_mu;        end
if isfield(opts, 'admm_tol');   admm_tol = opts.admm_tol;    end
if isfield(opts, 'admm_iter');	admm_iter = opts.admm_iter;  end
%U 设置
U = tensor_dc_px(X , theta);
err = 1;
for k = 1 : admm_iter%%%%%%%%%%%%%%%%主干有opts.admm_iter
    last_X = X;
    % update S
    [S, tnn, trank] = prox_capped_tnn(X + Y/mu, 1/mu/theta);
    %-----重要!!!------------------
    S(omega) = M(omega);
    % update X
    X = U/mu + S - Y/mu ;
    if mod(k, DISPLAY_EVERY) == 0
        errlast_X = norm(last_X(:));
        errX = norm(X(:));
        err = errX - errlast_X;        
        fprintf(['iter %d,\t mu=%.4f,\t t-rank=%d,' ...
            '\t obj=%.4f,\t err=%.4f\n'], k, mu, trank, tnn, err);
    end
    if (err < 0.1) && (err > -0.1)
        break;
    end
    % update Y
    Y = Y + mu * (X - S);
    % update mu
    mu = min( rho * mu, max_mu);
    % mu 有什么作用？？？
end

X_opt = X;
iter = k;

end