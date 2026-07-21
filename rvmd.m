function [mode, info, restart] = rvmd(Q, K, Alpha, varargin)
%RVMD Reduced-order variational mode decomposition.
%   [MODE, INFO] = RVMD(Q, K, ALPHA) decomposes the S-by-T data matrix Q
%   into K modes. Q may be real or complex. MODE.PHI is S-by-K,
%   MODE.C is T-by-K, and Q is approximated by MODE.PHI * MODE.C.'.
%
%   ALPHA is the bandwidth penalty in the RVMD filter
%       1 / (1 + 2*ALPHA*(|f| - f_k)^2).
%   Frequencies are reported in cycles per sample in [0, 0.5].
%
%   [...] = RVMD(Q, K, ALPHA, NAME, VALUE) accepts:
%       'Weight'           positive scalar or S-vector (default 1)
%       'Tolerance'        nonnegative stopping tolerance (default 5e-3)
%       'MaximumSteps'     positive total iteration cap (default 500)
%       'InitFreqType'     -1 random, 0 zero, 1 uniform (default 1)
%       'InitFreqMaximum'  initial upper frequency, capped at 0.5
%       'Device'           'cpu' or 'gpu' (default 'cpu')
%       'FPPrecision'      'single' or 'double' (default 'single')
%       'nDC'              number of leading modes fixed at zero frequency
%       'Display'          'off' or 'iter' (default 'off')
%       'Restart'          restart struct returned by an earlier call
%
%   [MODE, INFO, RESTART] = RVMD(...) also returns the internal state.
%   Resume to a total iteration cap N with
%       RVMD('Restart', RESTART, 'MaximumSteps', N).
%
%   Real Q uses the nonnegative half-spectrum and real spatial modes.
%   Complex Q uses the full spectrum, complex spatial modes, and the even
%   distance |f|-f_k. Complex modes are phase-normalized so the largest
%   entry of each spatial mode is real and positive.
%
%   Reference:
%   Liao et al. (2023), Journal of Fluid Mechanics 966, A7.
%   https://doi.org/10.1017/jfm.2023.435

if nargin >= 1 && isTextScalar(Q) && strcmpi(char(Q), 'Restart')
    if nargin < 2 || ~isstruct(K)
        error('rvmd:InvalidRestart', ...
            'RVMD(''Restart'', state, ...) requires a restart struct.');
    end
    restartArguments = varargin;
    if nargin >= 3
        restartArguments = [{Alpha}, varargin];
    end
    [settings, state] = parseRestartCall(K, restartArguments{:});
    Q = state.Q;
    K = state.K;
    Alpha = state.alpha;
    restarting = true;
else
    if nargin < 3
        error('rvmd:NotEnoughInputs', ...
            'Q, K, and Alpha are required for a new decomposition.');
    end
    [settings, restartInput] = parseNewCall(Q, K, Alpha, varargin{:});
    if isstruct(restartInput)
        restarting = true;
        [settings, state] = parseRestartCall(restartInput, ...
            'Tolerance', settings.Tolerance, ...
            'MaximumSteps', settings.MaximumSteps, ...
            'Device', settings.Device, 'Display', settings.Display);
        Q = state.Q;
        K = state.K;
        Alpha = state.alpha;
    elseif isequal(restartInput, 0)
        restarting = false;
        state = struct();
    else
        error('rvmd:InvalidRestart', ...
            'Restart must be zero or a restart struct returned by RVMD.');
    end
end

if restarting
    precision = state.FPPrecision;
    Q = cast(state.Q, precision);
    K = state.K;
    Alpha = cast(state.alpha, precision);
    weight = restoreWeight(state.weight, state.S, precision);
    isRealInput = state.isRealInput;
    nDC = state.nDC;
    phi = cast(state.phi_n, precision);
    coefficientSpectrum = cast(state.c_spec_n, precision);
    completedSteps = state.Iteration.steps;
    frequencyHistoryStored = cast(state.Iteration.omega, precision);
    differenceStored = cast(state.Iteration.difference, precision);
    initFreqType = state.InitFreqType;
    initFreqMaximum = state.InitFreqMaximum;
else
    precision = settings.FPPrecision;
    Q = cast(Q, precision);
    Alpha = cast(Alpha, precision);
    weight = normalizeWeight(settings.Weight, size(Q, 1), precision);
    isRealInput = isreal(Q);
    nDC = settings.nDC;
    completedSteps = 0;
    frequencyHistoryStored = zeros(K, 0, precision);
    differenceStored = zeros(1, 0, precision);
    initFreqType = settings.InitFreqType;
    initFreqMaximum = settings.InitFreqMaximum;
end

S = size(Q, 1);
T = size(Q, 2);
extendedLength = 2 * T;
half = ceil(T / 2);

QExtended = zeros(S, extendedLength, precision);
QExtended(:, 1:half) = Q(:, half:-1:1);
QExtended(:, (half + 1):(half + T)) = Q;
QExtended(:, (half + T + 1):end) = Q(:, T:-1:(half + 1));

QSpectrum = fft(QExtended, [], 2);
clear QExtended
if isRealInput
    spectrumLength = T + 1;
    QSpectrum = QSpectrum(:, 1:spectrumLength);
    frequency = cast((0:T).' / extendedLength, precision);
    frequencyWeight = ones(spectrumLength, 1, precision);
    frequencyWeight([1, end]) = cast(0.5, precision);
else
    spectrumLength = extendedLength;
    QSpectrum = fftshift(QSpectrum, 2);
    frequency = cast((-T:(T - 1)).' / extendedLength, precision);
    frequencyWeight = ones(spectrumLength, 1, precision);
end
absoluteFrequency = abs(frequency);

if restarting
    if size(phi, 1) ~= S || size(phi, 2) ~= K || ...
            size(coefficientSpectrum, 1) ~= spectrumLength || ...
            size(coefficientSpectrum, 2) ~= K
        error('rvmd:InvalidRestart', ...
            'Restart state dimensions do not match its saved data.');
    end
else
    phi = initialSpatialModes(S, K, weight, precision, isRealInput);
    if scalarValue(norm(Q, 'fro')) == 0
        coefficientSpectrum = zeros(spectrumLength, K, precision);
    else
        coefficientSpectrum = ones(spectrumLength, K, precision) ...
            * cast(eps(precision), precision);
    end
end

historyCapacity = max(settings.MaximumSteps + 1, completedSteps + 1);
frequencyHistory = zeros(K, historyCapacity, precision);
differenceHistory = zeros(1, max(settings.MaximumSteps, completedSteps), precision);
if restarting
    frequencyHistory(:, 1:(completedSteps + 1)) = ...
        frequencyHistoryStored(:, 1:(completedSteps + 1));
    if completedSteps > 0
        differenceHistory(1:completedSteps) = differenceStored(1:completedSteps);
    end
else
    frequencyHistory(:, 1) = initializeFrequencies( ...
        K, nDC, initFreqType, initFreqMaximum, precision);
end

if strcmp(settings.Device, 'gpu')
    if exist('gpuArray', 'file') ~= 2 && exist('gpuArray', 'class') ~= 8
        error('rvmd:GPUUnavailable', ...
            'GPU computation requires MATLAB Parallel Computing Toolbox.');
    end
    QSpectrum = gpuArray(QSpectrum);
    phi = gpuArray(phi);
    coefficientSpectrum = gpuArray(coefficientSpectrum);
    frequencyHistory = gpuArray(frequencyHistory);
    weight = gpuArray(weight);
    absoluteFrequency = gpuArray(absoluteFrequency);
    frequencyWeight = gpuArray(frequencyWeight);
end

if restarting && isfield(state, 'residual_n')
    residual = cast(state.residual_n, precision);
    if ~isequal(size(residual), size(QSpectrum))
        error('rvmd:InvalidRestart', ...
            'Saved residual dimensions do not match the restart state.');
    end
    if strcmp(settings.Device, 'gpu')
        residual = gpuArray(residual);
    end
else
    residual = QSpectrum - phi * coefficientSpectrum.';
end
dataSpectrumNorm = scalarValue(norm(QSpectrum, 'fro'));
clear QSpectrum
scaleFloor = cast(eps(precision) * max(dataSpectrumNorm, 1), precision);
if strcmp(settings.Device, 'gpu')
    scaleFloor = gpuArray(scaleFloor);
end
iteration = completedSteps + 1;
if completedSteps == 0
    difference = inf;
else
    difference = scalarValue(differenceHistory(completedSteps));
end

while iteration <= settings.MaximumSteps && difference > settings.Tolerance
    differenceAccumulator = zeros(1, 'like', frequencyHistory);

    for k = 1:K
        oldMode = phi(:, k) * coefficientSpectrum(:, k).';
        residual = residual + oldMode;

        projection = residual * ...
            (conj(coefficientSpectrum(:, k)) .* frequencyWeight);
        if isRealInput
            projection = real(projection);
        end
        projectionNorm = sqrt(sum(abs(projection).^2 .* weight));

        if isPositiveFinite(projectionNorm)
            phi(:, k) = projection / projectionNorm;
        else
            phi(:, k) = normalizedFallback( ...
                phi(:, k), k, weight, precision, isRealInput);
        end
        phi(:, k) = canonicalizePhase(phi(:, k), isRealInput);

        denominator = 1 + 2 * Alpha * ...
            (absoluteFrequency - frequencyHistory(k, iteration)).^2;
        coefficientSpectrum(:, k) = ...
            (residual.' * (conj(phi(:, k)) .* weight)) ./ denominator;

        coefficientEnergy = sum(frequencyWeight .* ...
            abs(coefficientSpectrum(:, k)).^2);
        if k <= nDC
            frequencyHistory(k, iteration + 1) = 0;
        elseif isPositiveFinite(coefficientEnergy)
            frequencyHistory(k, iteration + 1) = ...
                sum(frequencyWeight .* absoluteFrequency .* ...
                abs(coefficientSpectrum(:, k)).^2) / coefficientEnergy;
        else
            frequencyHistory(k, iteration + 1) = ...
                frequencyHistory(k, iteration);
        end

        newMode = phi(:, k) * coefficientSpectrum(:, k).';
        residual = residual - newMode;

        oldNorm = norm(oldMode, 'fro');
        changeNorm = norm(newMode - oldMode, 'fro');
        relativeChange = changeNorm / max(oldNorm, scaleFloor);
        differenceAccumulator = differenceAccumulator + relativeChange;
    end

    difference = scalarValue(differenceAccumulator);
    differenceHistory(iteration) = differenceAccumulator;
    if strcmp(settings.Display, 'iter')
        fprintf('iteration step: %d    difference: %.8g\n', ...
            iteration, difference);
    end
    iteration = iteration + 1;
end

completedSteps = iteration - 1;
if strcmp(settings.Device, 'gpu')
    phi = gather(phi);
    coefficientSpectrum = gather(coefficientSpectrum);
    frequencyHistory = gather(frequencyHistory);
    weight = gather(weight);
    residual = gather(residual);
end

if isRealInput
    fullSpectrum = zeros(extendedLength, K, 'like', coefficientSpectrum);
    fullSpectrum(1:spectrumLength, :) = coefficientSpectrum;
    fullSpectrum((spectrumLength + 1):end, :) = ...
        conj(coefficientSpectrum((spectrumLength - 1):-1:2, :));
    coefficient = real(ifft(fullSpectrum, [], 1));
    phi = real(phi);
else
    coefficient = ifft(ifftshift(coefficientSpectrum, 1), [], 1);
end
coefficient = coefficient((half + 1):(half + T), :);

finalFrequency = frequencyHistory(:, completedSteps + 1);
[frequencySorted, order] = sort(finalFrequency, 'ascend');
mode.phi = phi(:, order);
mode.c = coefficient(:, order);
mode.omega = frequencySorted;
mode.energy = sum(abs(mode.c).^2, 1);

info.S = S;
info.T = T;
info.K = K;
info.alpha = Alpha;
info.weight = weight;
info.Tolerance = settings.Tolerance;
info.MaximumSteps = settings.MaximumSteps;
info.InitFreqType = initFreqType;
info.InitFreqMaximum = initFreqMaximum;
info.Device = settings.Device;
info.FPPrecision = precision;
info.nDC = nDC;
info.isRealInput = isRealInput;
info.Iteration.steps = completedSteps;
info.Iteration.omega = frequencyHistory(:, 1:(completedSteps + 1));
info.Iteration.difference = differenceHistory(1:completedSteps);
if completedSteps == 0
    info.Iteration.converged = false;
else
    info.Iteration.converged = ...
        differenceHistory(completedSteps) <= settings.Tolerance;
end

restart.version = 3;
restart.Q = Q;
restart.S = S;
restart.T = T;
restart.K = K;
restart.alpha = Alpha;
restart.weight = weight;
restart.Tolerance = settings.Tolerance;
restart.MaximumSteps = settings.MaximumSteps;
restart.InitFreqType = initFreqType;
restart.InitFreqMaximum = initFreqMaximum;
restart.Device = settings.Device;
restart.FPPrecision = precision;
restart.nDC = nDC;
restart.isRealInput = isRealInput;
restart.Display = settings.Display;
restart.Iteration = info.Iteration;
restart.c_spec_n = coefficientSpectrum;
restart.phi_n = phi;
restart.residual_n = residual;
end

function [settings, restartInput] = parseNewCall(Q, K, Alpha, varargin)
p = inputParser;
p.FunctionName = 'rvmd';
p.PartialMatching = false;
addRequired(p, 'Q');
addRequired(p, 'K');
addRequired(p, 'Alpha');
addParameter(p, 'Restart', 0);
addParameter(p, 'Weight', 1);
addParameter(p, 'Tolerance', 5e-3);
addParameter(p, 'MaximumSteps', 500);
addParameter(p, 'InitFreqType', 1);
addParameter(p, 'InitFreqMaximum', 0.5);
addParameter(p, 'Device', 'cpu');
addParameter(p, 'FPPrecision', 'single');
addParameter(p, 'nDC', 0);
addParameter(p, 'Display', 'off');
parse(p, Q, K, Alpha, varargin{:});

validateProblem(Q, K, Alpha);
settings = validateSettings(p.Results, K);
restartInput = p.Results.Restart;
end

function [settings, state] = parseRestartCall(restartInput, varargin)
state = upgradeRestartState(restartInput);
p = inputParser;
p.FunctionName = 'rvmd';
p.PartialMatching = false;
addParameter(p, 'Tolerance', state.Tolerance);
addParameter(p, 'MaximumSteps', state.MaximumSteps);
addParameter(p, 'Device', 'cpu');
addParameter(p, 'Display', state.Display);
% Accepted for compatibility; restart state remains authoritative.
addParameter(p, 'InitFreqType', state.InitFreqType);
addParameter(p, 'InitFreqMaximum', state.InitFreqMaximum);
addParameter(p, 'FPPrecision', state.FPPrecision);
addParameter(p, 'nDC', state.nDC);
parse(p, varargin{:});

compatibility = p.Results;
if compatibility.InitFreqType ~= state.InitFreqType || ...
        compatibility.InitFreqMaximum ~= state.InitFreqMaximum || ...
        ~strcmpi(char(compatibility.FPPrecision), state.FPPrecision) || ...
        compatibility.nDC ~= state.nDC
    error('rvmd:ImmutableRestartOption', ...
        'Initialization, precision, and nDC cannot change on restart.');
end

settings = struct();
settings.Tolerance = validateNonnegativeScalar( ...
    p.Results.Tolerance, 'Tolerance');
settings.MaximumSteps = validatePositiveInteger( ...
    p.Results.MaximumSteps, 'MaximumSteps');
if settings.MaximumSteps < state.Iteration.steps
    error('rvmd:InvalidMaximumSteps', ...
        'MaximumSteps cannot be less than the completed restart steps.');
end
settings.InitFreqType = state.InitFreqType;
settings.InitFreqMaximum = state.InitFreqMaximum;
settings.Device = validateChoice(p.Results.Device, ...
    {'cpu', 'gpu'}, 'rvmd:InvalidDevice', 'Device');
settings.FPPrecision = state.FPPrecision;
settings.nDC = state.nDC;
settings.Display = validateChoice(p.Results.Display, ...
    {'off', 'iter'}, 'rvmd:InvalidDisplay', 'Display');
end

function validateProblem(Q, K, Alpha)
if ~isnumeric(Q) || isempty(Q) || ndims(Q) > 2 || ...
        ~all(isfinite(Q(:)))
    error('rvmd:InvalidData', ...
        'Q must be a nonempty, finite numeric matrix.');
end
validatePositiveInteger(K, 'K');
if ~isnumeric(Alpha) || ~isreal(Alpha) || ~isscalar(Alpha) || ...
        ~isfinite(Alpha) || Alpha < 0
    error('rvmd:InvalidAlpha', ...
        'Alpha must be a finite, nonnegative real scalar.');
end
end

function settings = validateSettings(results, K)
settings = struct();
settings.Weight = results.Weight;
settings.Tolerance = validateNonnegativeScalar( ...
    results.Tolerance, 'Tolerance');
settings.MaximumSteps = validatePositiveInteger( ...
    results.MaximumSteps, 'MaximumSteps');
if ~isnumeric(results.InitFreqType) || ~isscalar(results.InitFreqType) || ...
        ~ismember(results.InitFreqType, [-1, 0, 1])
    error('rvmd:InvalidInitFreqType', ...
        'InitFreqType must be -1, 0, or 1.');
end
settings.InitFreqType = double(results.InitFreqType);
settings.InitFreqMaximum = min(validateNonnegativeScalar( ...
    results.InitFreqMaximum, 'InitFreqMaximum'), 0.5);
settings.Device = validateChoice(results.Device, {'cpu', 'gpu'}, ...
    'rvmd:InvalidDevice', 'Device');
settings.FPPrecision = validateChoice(results.FPPrecision, ...
    {'single', 'double'}, 'rvmd:InvalidPrecision', 'FPPrecision');
settings.nDC = validateNonnegativeInteger(results.nDC, 'nDC');
if settings.nDC > K
    error('rvmd:InvalidNDC', 'nDC must be no greater than K.');
end
settings.Display = validateChoice(results.Display, {'off', 'iter'}, ...
    'rvmd:InvalidDisplay', 'Display');
end

function state = upgradeRestartState(input)
required = {'Q', 'S', 'T', 'K', 'alpha', 'weight', ...
    'c_spec_n', 'phi_n'};
if ~isstruct(input) || ~all(isfield(input, required))
    error('rvmd:InvalidRestart', ...
        'The restart struct is missing required RVMD state fields.');
end
state = input;

if isfield(input, 'Iteration')
    state.Iteration = input.Iteration;
elseif all(isfield(input, {'steps', 'omega', 'difference'}))
    state.Iteration.steps = input.steps;
    state.Iteration.omega = input.omega;
    state.Iteration.difference = input.difference;
else
    error('rvmd:InvalidRestart', ...
        'The restart struct has no complete iteration history.');
end

if ~isfield(state, 'FPPrecision')
    if isa(state.Q, 'single')
        state.FPPrecision = 'single';
    else
        state.FPPrecision = 'double';
    end
else
    state.FPPrecision = validateChoice(state.FPPrecision, ...
        {'single', 'double'}, 'rvmd:InvalidRestart', 'FPPrecision');
end
if ~isfield(state, 'isRealInput')
    state.isRealInput = size(state.c_spec_n, 1) == state.T + 1;
end
if ~isfield(state, 'nDC')
    state.nDC = 0;
end
if ~isfield(state, 'Tolerance')
    state.Tolerance = 5e-3;
end
if ~isfield(state, 'MaximumSteps')
    state.MaximumSteps = max(500, state.Iteration.steps);
end
if ~isfield(state, 'InitFreqType')
    state.InitFreqType = 1;
end
if ~isfield(state, 'InitFreqMaximum')
    state.InitFreqMaximum = 0.5;
end
if ~isfield(state, 'Display')
    state.Display = 'off';
end

if state.S ~= size(state.Q, 1) || state.T ~= size(state.Q, 2) || ...
        state.K ~= size(state.phi_n, 2) || ...
        state.Iteration.steps < 0 || ...
        size(state.Iteration.omega, 1) ~= state.K || ...
        size(state.Iteration.omega, 2) < state.Iteration.steps + 1 || ...
        numel(state.Iteration.difference) < state.Iteration.steps
    error('rvmd:InvalidRestart', 'The restart struct is internally inconsistent.');
end
validateProblem(state.Q, state.K, state.alpha);
state.nDC = validateNonnegativeInteger(state.nDC, 'nDC');
if state.nDC > state.K
    error('rvmd:InvalidRestart', 'Saved nDC is greater than saved K.');
end
state.Tolerance = validateNonnegativeScalar(state.Tolerance, 'Tolerance');
state.MaximumSteps = validatePositiveInteger( ...
    state.MaximumSteps, 'MaximumSteps');
state.Display = validateChoice(state.Display, {'off', 'iter'}, ...
    'rvmd:InvalidRestart', 'Display');
end

function weight = normalizeWeight(inputWeight, S, precision)
weight = expandWeight(inputWeight, S);
weight = cast(weight / mean(weight), precision);
end

function weight = restoreWeight(inputWeight, S, precision)
weight = cast(expandWeight(inputWeight, S), precision);
end

function weight = expandWeight(inputWeight, S)
if ~isnumeric(inputWeight) || ~isreal(inputWeight) || ...
        isempty(inputWeight) || ~all(isfinite(inputWeight(:))) || ...
        any(inputWeight(:) <= 0) || ...
        ~(isscalar(inputWeight) || numel(inputWeight) == S)
    error('rvmd:InvalidWeight', ...
        'Weight must be a positive finite scalar or a vector of length S.');
end
if isscalar(inputWeight)
    weight = repmat(inputWeight, S, 1);
else
    weight = inputWeight(:);
end
end

function phi = initialSpatialModes(S, K, weight, precision, isRealInput)
phi = zeros(S, K, precision);
for k = 1:K
    index = mod(k - 1, S) + 1;
    phi(index, k) = 1 / sqrt(weight(index));
end
if ~isRealInput
    phi = complex(phi);
end
end

function omega = initializeFrequencies( ...
        K, nDC, initFreqType, initFreqMaximum, precision)
omega = zeros(K, 1, precision);
freeModes = K - nDC;
if freeModes == 0
    return
end
switch initFreqType
    case -1
        omega((nDC + 1):end) = cast( ...
            rand(freeModes, 1) * initFreqMaximum, precision);
    case 0
        % Already zero.
    case 1
        omega((nDC + 1):end) = cast( ...
            (1:freeModes).' / freeModes * initFreqMaximum, precision);
end
end

function phi = normalizedFallback(phi, k, weight, precision, isRealInput)
phiNorm = sqrt(sum(abs(phi).^2 .* weight));
if isPositiveFinite(phiNorm)
    phi = phi / phiNorm;
else
    phi = zeros(size(phi), 'like', phi);
    index = mod(k - 1, numel(phi)) + 1;
    phi(index) = 1 / sqrt(weight(index));
end
if isRealInput
    phi = real(phi);
end
end

function phi = canonicalizePhase(phi, isRealInput)
[~, pivot] = max(abs(phi));
pivot = scalarValue(pivot);
pivotValue = phi(pivot);
if isRealInput
    if scalarValue(pivotValue) < 0
        phi = -phi;
    end
else
    pivotMagnitude = abs(pivotValue);
    if isPositiveFinite(pivotMagnitude)
        phi = phi * conj(pivotValue) / pivotMagnitude;
        phi(pivot) = complex(pivotMagnitude, 0);
    end
end
end

function value = validateNonnegativeScalar(value, name)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 0
    error(['rvmd:Invalid', name], ...
        '%s must be a finite nonnegative real scalar.', name);
end
value = double(value);
end

function value = validatePositiveInteger(value, name)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 1 || value ~= fix(value)
    error(['rvmd:Invalid', name], '%s must be a positive integer.', name);
end
value = double(value);
end

function value = validateNonnegativeInteger(value, name)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 0 || value ~= fix(value)
    error(['rvmd:Invalid', upper(name)], ...
        '%s must be a nonnegative integer.', name);
end
value = double(value);
end

function value = validateChoice(value, choices, identifier, name)
if ~isTextScalar(value)
    error(identifier, '%s must be a text scalar.', name);
end
value = lower(char(value));
if ~any(strcmp(value, choices))
    error(identifier, '%s has an unsupported value.', name);
end
end

function tf = isTextScalar(value)
tf = ischar(value) && (isrow(value) || isempty(value));
if ~tf && (exist('isstring', 'builtin') == 5 || ...
        exist('isstring', 'file') == 2)
    tf = isstring(value) && isscalar(value);
end
end

function tf = isPositiveFinite(value)
value = scalarValue(value);
tf = isfinite(value) && value > 0;
end

function value = scalarValue(value)
if isa(value, 'gpuArray')
    value = gather(value);
end
value = double(value);
end
