function [phi, c, omega, energy, omega_iter, diff_iter] = rvmd(Q, K, alpha, varargin)
% Reduced-order variational mode decomposition (RVMD)
% Description
%   [phi,c,omega,energy,omega_iter,diff_iter] = rvmd(Q,K,alpha) returns the
%   reduced-order variational mode decomposition results of the data matrix
%   Q whose first dimension is space, i.e., size(Q)=[S, T], where S and T
%   are the numbers of sampling points in space and time, respectively.
%   The columns of phi are the spatial distributions, and the columns of c
%   are the corresponding time-evolution coefficients. Vector omega and
%   energy are the central frequencies and the mode energy, respectively.
%   omega_iter and diff_iter contain the iteration values of the central
%   frequencies and the iteration difference. All the components of RVMD
%   modes are real-valued. At least two inputs are required: the number of
%   RVMD modes K, the filtering parameter alpha.
%
%   [...] = rvmd(Q,K,alpha,tol,N) specifies the termination criterions. The
%   loop will stop when the iteration difference is less than the given
%   tolerance tol or when the number of iteration steps reaches N.
%
%   [...] = rvmd(Q,K,alpha,tol,N,init,initFreqMax) specifies how to
%   initialize the central frequencies in domain [0, min(initFreqMax, 0.5)]:
%       init            Strategy
%       -1              randomly distributed
%       0               all zero
%       1               uniformly distributed
%
%   Default values:
%       tol             5e-3
%       N               500
%       init            1
%       initFreqMax     0.5
% Reference
%   [1] Liao, Z.-M., Zhao, Z., Chen, L.-B., Wan, Z.-H., Liu, N.-S.
%       & Lu, X.-Y. 2023 Reduced-order variational mode decomposition:
%       reveal transient and nonstationary dynamics in fluid flows.
%       J. Fluid Mech. In press.


%% pre-processing
% mirror extension
Q_origin = Q; Q = [];
S = size(Q_origin, 1); % spatial dimension
T_origin = size(Q_origin, 2); % original temporal dimension

Q(:, 1:ceil(T_origin / 2)) = Q_origin(:, ceil(T_origin / 2):-1:1);
Q(:, ceil(T_origin / 2 + 1):3 * ceil(T_origin / 2)) = Q_origin;
Q(:, ceil(3 * T_origin / 2 + 1):2 * T_origin) = ...
    Q_origin(:, T_origin:-1:ceil(T_origin / 2 + 1));

T = size(Q, 2); % temporal dimension

% parameters initialization
tol = []; N = []; init = []; initFreqMax = [];
numVarargin = length(varargin);

if numVarargin >= 1
    tol = varargin{1};
    
    if numVarargin >= 2
        N = varargin{2};
        
        if numVarargin >= 3
            init = varargin{3};
            
            if numVarargin >= 4
                initFreqMax = min(varargin{4}, 0.5);
            end
            
        end
        
    end
    
end

if isempty(initFreqMax) % set default values
    initFreqMax = 0.5;
    
    if isempty(init)
        init = 1;
        
        if isempty(N)
            N = 500;
            
            if isempty(tol)
                tol = 5e-3;
            end
            
        end
        
    end
    
end

% variables initialization
Q_spec = fft(Q, [], 2);
T_spec = floor(T / 2 + 1);
Q_spec = Q_spec(:, 1:T_spec); % positive half in the frequency domain

phi_n = zeros(S, K) + eps; % spatial distributions
c_spec_n = zeros(T_spec, K) + eps; % time-evolution coefficients
mode_k = zeros(S, T_spec);
residual_k = zeros(S, T_spec);
omega_k = zeros(K, N); % central frequencies

switch init
    case - 1
        omega_k(:, 1) = rand(K, 1) * initFreqMax;
    case 0
        omega_k(:, 1) = 0;
    case 1
        omega_k(:, 1) = (0:1 / (K - 1):1) * initFreqMax;
end

omega = (0:floor(T / 2)) / T; % frequency list

%% main loop
diff = tol + eps;

n = 1;
residual_k = Q_spec - phi_n * c_spec_n.';

while (n <= N && diff > tol)
    
    diff = 0;
    
    for k = 1:K
        % calculate residual matrix
        mode_k = phi_n(:, k) * c_spec_n(:, k).';
        residual_k = residual_k + mode_k;
        
        % update phi_k
        phi_n(:, k) = real(conj(residual_k) * c_spec_n(:, k));
        phi_n(:, k) = phi_n(:, k) / norm(phi_n(:, k), 'fro');
        
        % update c_k
        c_spec_n(:, k) = (residual_k.' * phi_n(:, k)) ./ ...
            (1 + 2 * alpha * (omega - omega_k(k, n)).^2).';
        
        % update omega_k
        omega_k(k, n + 1) = omega * (conj(c_spec_n(:, k)) .* c_spec_n(:, k)) ...
            / norm(c_spec_n(:, k), 'fro')^2;
        
        mode_k_new = phi_n(:, k) * c_spec_n(:, k).';
        residual_k = residual_k - mode_k_new;
        
        % calculate iteration difference
        diff = diff + norm(mode_k_new - mode_k, 'fro') / norm(mode_k, 'fro');
    end
    
    diff_iter(n) = diff;
    disp(['iteration step: ', num2str(n), '    differences: ', num2str(diff)]);
    
    n = n + 1;
end

%% post-processing
% central frequencies
omega_iter = omega_k(:, 1:n);
omega = omega_k(:, n);

% reconstruct time coefficients
c_spec = c_spec_n;
c_spec(T_spec + 1:T, :) = conj(c_spec_n(end - 1:-1:2, :));
c = ifft(c_spec);
c = c(T / 4 + 1:3 * T / 4, :);
% c_spec=fft(c,[],1);

% sort the RVMD modes from low frequency to high frequency
[~, index] = sort(omega);

for k = 1:K
    phi(:, k) = phi_n(:, index(k));
    c_new(:, k) = c(:, index(k));
    %     c_spec_new(:,index(k))=c_spec(:,index(k));
    omega_new(k, :) = omega(index(k), :);
    energy(k) = norm(c(:, index(k)))^2;
end

c = c_new;
% c_spec=c_spec_new;
omega = omega_new;

end