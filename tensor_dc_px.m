function p_x = tensor_dc_px(X , theta)

% the singular value decomposition of a 3D tensor

dim = ndims(X);
[n1, n2, n3] = size(X);
n12 = min(n1, n2);
Xf = fft(X, [], dim);
Uf = zeros(n1, n12, n3);
Vf = zeros(n2, n12, n3);
Sf = zeros(n12, n12, n3);%
%%%
Bf = Sf;%
%%%
Xf(isnan(Xf)) = 0;
Xf(isinf(Xf)) = 0;
%
%theta = 1;% theta取值值得商榷！cappedtnnadmm
%
t_rank = 0;
for i = 1 : n3
    [Uf(:,:,i), Sf(:,:,i), Vf(:,:,i)] = svd(Xf(:,:,i), 'econ');
    diagS = diag(Sf(:, :, i));%向量化了的diagS --- a vector
    diagB = diagS;%向量化了的diagB --- a vector
    for j = 1 : length(diagB)
        if diagB(j) > theta
            diagB(j) = 1;
        elseif diagB(i) < theta
            diagB(j) = 0;
        else
            diagB(j) = 0.5;
        end
    end
    Bf(:, :, i) = diag(diagB);%矩阵化了的diagB
    % if xxx > theta b(i,i) = 1 end
    % else if xxx < theta b(i,i) = 0 end
    % else b(i,i) = 0.5
    % 针对S进行改造
    temp = length(find(diagS > 0));
    t_rank = max(temp, t_rank);
end

Uf = Uf(:, 1:t_rank, :);
%Sf = Sf(1:t_rank, 1:t_rank, :);
Vf = Vf(:, 1:t_rank, :);
Bf = Bf(1:t_rank, 1:t_rank, :);


U = ifft(Uf, [], dim);
%S = ifft(Sf, [], dim);
V = ifft(Vf, [], dim);
B = ifft(Bf, [], dim);

%temp1 = tprod(U, B);
p_x = tprod( tprod( U, B), tran(V) );
end