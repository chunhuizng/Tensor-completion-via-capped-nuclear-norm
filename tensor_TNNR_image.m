%% add path
addpath(genpath(cd))
close all
clear
clc

%% read image files directory information
admm_result = './result/admm/image';
apgl_result = './result/apgl/image';
if ~exist(admm_result, 'dir'),	mkdir(admm_result);	end
if ~exist(apgl_result, 'dir'),  mkdir(apgl_result);	end
% image_list = {'re1.jpg', 're2.jpg', 're3.jpg', 're4.jpg', 're5.jpg', ...
%               're6.jpg', 're7.jpg', 're8.jpg', 're9.jpg', 're10.jpg', ...
%              };
image_list = {'new1.jpg', 'new2.jpg', 'new3.jpg', 'new4.jpg', 'new5.jpg', ...
              'new6.jpg', 'new7.jpg', 'new8.jpg', 'new9.jpg', 'new10.jpg', ...
             };
         
file_list = dir('mask');
num_mask = length(file_list) - 2;
mask_list = cell(num_mask+2, 1);
for i = 1 : num_mask
    mask_list{i} = file_list(i+2).name; 
end

%% parameter configuration
image_id = 1;           % select an image for experiment
mask_id  = 4;           % select a mask for experiment

opts.block = 0;         % 1 for block occlusion, 0 for random noise
opts.lost = 0.50;       % percentage of lost elements in matrix
opts.save_eps = 1;      % save eps figure in result directory
% it requires to test all ranks from min_R to max_R, note that different
% images have different ranks, and various masks affect the ranks, too.
%^^^^^^^^^^^^^^^^^^^^^^^^^
opts.min_R = 1 ;       %控制theta值范围 theta = 0.01 * R ;
opts.max_R = 1;        %控制theta值范围 theta = 0.01 * R ;
%__________________________
opts.out_iter = 50;   % !!!!!!!!!!!!!!!!! outer iteration 50
opts.out_tol = 1e-3;    % tolerance of outer iteration

%什么意思？
opts.mu = 6.5e-3;         % mu of ADMM optimization
opts.rho = 1.25;        % rho of ADMM optimization 1.05
opts.max_mu = 1e10;     % max value of mu 
opts.admm_tol = 1e-4;   % tolerance of ADMM iteration
%什么意思？
opts.admm_iter = 200;   % !!!!!!!!!!!!!!!!!ADMM iteration 200
opts.maxP = 255;        % max pixel value

%% select an image and a mask for experiment
image_name = image_list{image_id};
X_full = double(imread(image_name));
[n1, n2, n3] = size(X_full);
fprintf('choose image: %s, ', image_name);

if opts.block  
    % block occlusion
    mask = double(imread(mask_list{mask_id}));
    mask = mask ./ max(mask(:));       % index matrix of the known elements
    fprintf('mask: %s.\n', mask_list{mask_id});
    omega = find(mask);% 定义Ω，mask是
else
    lost = opts.lost;
    fprintf('loss: %d%% elements are randomly missing\n', lost*100);

    mask = double(rand(n1,n2,n3) < (1-lost));
    omega = find(mask);
end

M = zeros(n1, n2, n3);
M(omega) = X_full(omega);% M is the incompelte image.
max_P = opts.maxP;       % maxP is the max pixel value.

%% tensor truncated tensor nuclear norm, using ADMM
fprintf('ADMM method to recover an image with missing pixels\n');
opts.method = 'ADMM';

t1 = tic;
[X_hat, admm_res] = capped_tensor_tnnr(X_full, omega, opts);
toc(t1)

admm_R = admm_res.best_R;
admm_psnr = admm_res.best_psnr;
admm_erec = admm_res.best_erec;
admm_time_cost = admm_res.time(admm_R);
admm_iteration = admm_res.iterations(admm_R);
admm_total_iter = admm_res.total_iter(admm_R);

figure
subplot(1,3,1)
imshow(X_full/max_P)
title('original image')
subplot(1,3,2)
imshow(M/max_P)
title('incomplete image')
subplot(1,3,3)
imshow(X_hat/max_P)
title('recovered image')

%% save eps figure in result directory
if opts.save_eps
    fig_eps = figure;
    imshow(X_hat ./ 255, 'border', 'tight');
    split_name = regexp(image_name, '[.]', 'split');
    fig_name = sprintf('%s/%s_rank_%d_PSNR_%.2f', ...
        admm_result, split_name{1}, admm_R, admm_psnr);
    saveas(gcf, [fig_name '.eps'], 'psc2');
    fprintf('eps figure saved in %s.eps\n', fig_name);
    close(fig_eps);
end

fprintf('\nCapped Tensor TNNR (ADMM):\n');
fprintf('theta=0.0%d, psnr=%.4f, erec=%.4f, time=%.3f s, iteration=%d(%d)\n', ...
    admm_R, admm_psnr, admm_erec, admm_time_cost, admm_iteration, ...
    admm_total_iter);

disp(' ');

figure('NumberTitle', 'off', 'Name', 'Capped Tensor TNNR (ADMM) result')
subplot(2, 2, 1)
plot(admm_res.R, admm_res.Psnr, 'o-')%%%要改admm_res.Psnr
xlabel('theta*10^-2')
ylabel('PSNR')

subplot(2, 2, 2)
plot(admm_res.R * 0.1, admm_res.Erec, 'diamond-')%%%要改admm_res.Erec
xlabel('theta*10^-2')
ylabel('Recovery error')

subplot(2, 2, 3)
plot(admm_res.Psnr_iter, 'square-')
xlabel('Iteration')
ylabel('PSNR')

subplot(2, 2, 4)
plot(admm_res.Erec_iter, '^-')
xlabel('Iteration')
ylabel('Recovery error')

%% record test results
outputFileName = fullfile(admm_result, 'parameters.txt'); 
fid = fopen(outputFileName, 'a') ;
fprintf(fid, '****** %s ******\n', datestr(now,0));
fprintf(fid, '%s\n', ['image: '           image_name               ]);
%fprintf(fid, '%s\n', ['mask: '            mask_list{mask_id}       ]);
fprintf(fid, '%s\n', ['block or noise: '  num2str(opts.block)      ]);
fprintf(fid, '%s\n', ['loss ratio: '      num2str(opts.lost)       ]);
fprintf(fid, '%s\n', ['save eps figure: ' num2str(opts.save_eps)   ]);
fprintf(fid, '%s\n', ['min rank: '        num2str(opts.min_R)      ]);
fprintf(fid, '%s\n', ['max rank: '        num2str(opts.max_R)      ]);
fprintf(fid, '%s\n', ['max iteration: '   num2str(opts.out_iter)   ]);
fprintf(fid, '%s\n', ['tolerance: '       num2str(opts.out_tol)    ]);
fprintf(fid, '%s\n', ['ADMM mu: '         num2str(opts.mu)         ]);
fprintf(fid, '%s\n', ['ADMM rho: '        num2str(opts.rho)        ]);
fprintf(fid, '%s\n', ['ADMM max_mu: '     num2str(opts.max_mu)     ]);
fprintf(fid, '%s\n', ['ADMM iteration: '  num2str(opts.admm_iter)  ]);
fprintf(fid, '%s\n', ['ADMM tolerance: '  num2str(opts.admm_tol)   ]);
fprintf(fid, '%s\n', ['max pixel value: ' num2str(opts.maxP)       ]);

fprintf(fid, '%s\n', ['rank: '            num2str(admm_R)       ]);
fprintf(fid, '%s\n', ['psnr: '            num2str(admm_psnr)       ]);
fprintf(fid, '%s\n', ['recovery error: '  num2str(admm_erec)       ]);
fprintf(fid, '%s\n', ['time cost: '       num2str(admm_time_cost)  ]);
fprintf(fid, 'iteration: %d(%d)\n',       admm_iteration, admm_total_iter);
fprintf(fid, '--------------------\n');
fclose(fid);