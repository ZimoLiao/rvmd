function [mode, info, restart] = rvmd(Q, K, Alpha, varargin)
%RVMD Reduced-order variational mode decomposition.
%   [MODE, INFO] = RVMD(Q, K, ALPHA) decomposes the finite S-by-T numeric
%   matrix Q into K modes. Q may be real or complex, K is a positive
%   integer, and ALPHA is a finite nonnegative scalar. The reconstruction
%   convention is
%       Q_APPROX = MODE.PHI * MODE.C.'
%   where .' is the nonconjugate transpose.
%
%   ALPHA is the bandwidth penalty in the RVMD filter
%       1 / (1 + 2*ALPHA*(|f| - f_k)^2).
%   Center frequencies are reported in cycles per sample in [0, 0.5].
%
%   [...] = RVMD(Q, K, ALPHA, NAME, VALUE) accepts:
%       'Weight'             Positive scalar or S-vector. It is normalized
%                            to have mean 1. Default: 1.
%       'Tolerance'          Finite nonnegative stopping tolerance applied
%                            to the sum of per-mode relative Frobenius
%                            changes over one sweep. Default: 5e-3.
%       'MaximumSteps'       Positive integer total sweep limit, including
%                            completed restart sweeps. Default: 500.
%       'InitFreqType'       For the M=K-nDC free modes: -1 samples uniform
%                            random values below InitFreqMaximum, 0 uses
%                            zeros, and 1 uses (1:M)/M*InitFreqMaximum.
%                            Ignored when InitialFrequencies is supplied.
%                            Default: 1.
%       'InitFreqMaximum'    Finite nonnegative initialization upper bound,
%                            capped at 0.5. Default: 0.5.
%       'Device'             'cpu' or 'gpu'. GPU execution requires MATLAB
%                            gpuArray support. Default: 'cpu'.
%       'FPPrecision'        'single' or 'double'. Default: 'single'.
%       'nDC'                Integer in [0,K]. The first nDC internal mode
%                            frequencies remain zero. Default: 0.
%       'InitialFrequencies' Empty or a finite K-vector in [0,0.5]. Its
%                            first nDC entries must be zero. Default: [].
%       'Display'            'off' prints nothing, 'final' prints the
%                            termination message, and 'iter' also prints
%                            sweep difference and sorted center frequencies.
%                            Default: 'off'.
%       'DisplayInterval'    Positive integer interval used by 'iter'
%                            display. Default: 20.
%       'OutputFcn'          Empty or STOP = FCN(PROGRESS,PHASE). Default: [].
%       'TimeLimit'          Nonnegative per-call seconds or Inf, checked
%                            before each complete sweep. Default: Inf.
%       'CheckpointFile'     MAT-file path. Empty disables automatic
%                            checkpointing. Default: ''.
%       'CheckpointInterval' Positive integer interval in completed sweeps.
%                            Default: 50.
%
%   MODE contains:
%       phi       S-by-K spatial modes, sorted by center frequency.
%       c         T-by-K time coefficients in the same order.
%       omega     K-by-1 sorted center frequencies, cycles per sample.
%       energy    1-by-K values SUM(ABS(c).^2,1).
%
%   INFO records the effective settings and termination status. ExitFlag is
%   1 for convergence, 0 for MaximumSteps, -1 for OutputFcn, and -2 for
%   TimeLimit. StopReason is respectively 'converged', 'maximumSteps',
%   'outputFunction', or 'timeLimit'. INFO.Iteration contains steps, the
%   unsorted K-by-(steps+1) omega history, the 1-by-steps difference
%   history, and converged.
%
%   RESTART is an opaque state for continuation. Resume to a total sweep
%   limit N with
%       RVMD('Restart', RESTART, 'MaximumSteps', N).
%   On restart, Tolerance, MaximumSteps, Device, Display, DisplayInterval,
%   OutputFcn, TimeLimit, CheckpointFile, and CheckpointInterval may be
%   changed. Initialization, precision, nDC, weights, K, and ALPHA remain
%   those of the saved state.
%
%   OutputFcn phases are 'init', 'iter', and 'done'. PROGRESS contains step,
%   difference, unsorted omega, elapsedTime, and stopReason. The callback
%   must return [] or a finite logical or numeric scalar. A nonzero value
%   during 'init' or 'iter' requests a stop at a consistent sweep boundary;
%   the value returned during 'done' is validated but not otherwise used.
%
%   CheckpointFile is written periodically and the final completed step is
%   also saved after normal termination. The MAT-file variable is named
%   restart. When a checkpoint already exists, it is retained at
%   [CheckpointFile '.prev'] before the new file is activated.
%
%   Real Q uses the nonnegative half-spectrum and real spatial modes.
%   Complex Q uses the full spectrum, complex spatial modes, and the even
%   distance |f|-f_k. Complex modes are phase-normalized so the largest
%   entry of each spatial mode is real and positive.
%
%   Reference:
%   Liao et al. (2023), Journal of Fluid Mechanics 966, A7.
%   https://doi.org/10.1017/jfm.2023.435

callTimer = tic;

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
        error('rvmd:RestartSyntax', ...
            ['Use RVMD(''Restart'', state, ...) to continue a saved ', ...
             'decomposition. Q, K, and Alpha are not accepted on restart.']);
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
    initialFrequencies = cast(state.InitialFrequencies, precision);
    S = state.S;
    T = state.T;
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
    S = size(Q, 1);
    T = size(Q, 2);
    if isempty(settings.InitialFrequencies)
        initialFrequencies = initializeFrequencies( ...
            K, nDC, initFreqType, initFreqMaximum, precision);
    else
        initialFrequencies = cast(settings.InitialFrequencies, precision);
    end
end

extendedLength = 2 * T;
half = ceil(T / 2);
if isRealInput
    spectrumLength = T + 1;
    frequency = cast((0:T).' / extendedLength, precision);
    frequencyWeight = ones(spectrumLength, 1, precision);
    frequencyWeight([1, end]) = cast(0.5, precision);
else
    spectrumLength = extendedLength;
    frequency = cast((-T:(T - 1)).' / extendedLength, precision);
    frequencyWeight = ones(spectrumLength, 1, precision);
end
absoluteFrequency = abs(frequency);

hasFastRestart = restarting && isfield(state, 'residual_n') && ...
    isfield(state, 'dataSpectrumNorm');
if hasFastRestart
    QSpectrum = [];
    dataSpectrumNorm = double(state.dataSpectrumNorm);
else
    if restarting
        if ~isfield(state, 'Q')
            error('rvmd:InvalidRestart', ...
                'This legacy restart state requires its saved Q field.');
        end
        sourceData = cast(state.Q, precision);
    else
        sourceData = Q;
    end
    QExtended = zeros(S, extendedLength, precision);
    QExtended(:, 1:half) = sourceData(:, half:-1:1);
    QExtended(:, (half + 1):(half + T)) = sourceData;
    QExtended(:, (half + T + 1):end) = ...
        sourceData(:, T:-1:(half + 1));
    QSpectrum = fft(QExtended, [], 2);
    clear QExtended sourceData
    if isRealInput
        QSpectrum = QSpectrum(:, 1:spectrumLength);
    else
        QSpectrum = fftshift(QSpectrum, 2);
    end
    dataSpectrumNorm = scalarValue(norm(QSpectrum, 'fro'));
end

if restarting
    if size(phi, 1) ~= S || size(phi, 2) ~= K || ...
            size(coefficientSpectrum, 1) ~= spectrumLength || ...
            size(coefficientSpectrum, 2) ~= K
        error('rvmd:InvalidRestart', ...
            'Restart state dimensions do not match its saved data.');
    end
else
    phi = initialSpatialModes(S, K, weight, precision, isRealInput);
    if dataSpectrumNorm == 0
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
    frequencyHistory(:, 1) = initialFrequencies;
end

if restarting && isfield(state, 'residual_n')
    residual = cast(state.residual_n, precision);
    if ~isequal(size(residual), [S, spectrumLength])
        error('rvmd:InvalidRestart', ...
            'Saved residual dimensions do not match the restart state.');
    end
else
    residual = QSpectrum - phi * coefficientSpectrum.';
end
clear QSpectrum Q

if strcmp(settings.Device, 'gpu')
    if exist('gpuArray', 'file') ~= 2 && exist('gpuArray', 'class') ~= 8
        error('rvmd:GPUUnavailable', ...
            'GPU computation requires MATLAB Parallel Computing Toolbox.');
    end
    phi = gpuArray(phi);
    coefficientSpectrum = gpuArray(coefficientSpectrum);
    frequencyHistory = gpuArray(frequencyHistory);
    weight = gpuArray(weight);
    absoluteFrequency = gpuArray(absoluteFrequency);
    frequencyWeight = gpuArray(frequencyWeight);
    residual = gpuArray(residual);
end
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

outputStopRequested = false;
timeLimitReached = false;
lastDisplayedStep = -1;
lastCheckpointStep = -1;
if ~isempty(settings.OutputFcn)
    progress = makeProgress(completedSteps, difference, ...
        frequencyHistory(:, completedSteps + 1), toc(callTimer));
    outputStopRequested = invokeOutputFcn( ...
        settings.OutputFcn, progress, 'init');
end

while iteration <= settings.MaximumSteps && ...
        difference > settings.Tolerance && ~outputStopRequested
    if toc(callTimer) >= settings.TimeLimit
        timeLimitReached = true;
        break
    end
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
                phi(:, k), k, weight, isRealInput);
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
    if strcmp(settings.Display, 'iter') && ...
            (iteration == completedSteps + 1 || ...
             mod(iteration, settings.DisplayInterval) == 0)
        displayIteration(iteration, difference, ...
            frequencyHistory(:, iteration + 1), lastDisplayedStep < 0);
        lastDisplayedStep = iteration;
    end
    if ~isempty(settings.OutputFcn)
        progress = makeProgress(iteration, difference, ...
            frequencyHistory(:, iteration + 1), toc(callTimer));
        outputStopRequested = invokeOutputFcn( ...
            settings.OutputFcn, progress, 'iter');
    end
    if ~isempty(settings.CheckpointFile) && ...
            mod(iteration, settings.CheckpointInterval) == 0
        checkpoint = makeRestartState(settings, S, T, K, Alpha, weight, ...
            isRealInput, precision, nDC, initFreqType, initFreqMaximum, ...
            initialFrequencies, iteration, frequencyHistory, ...
            differenceHistory, coefficientSpectrum, phi, residual, ...
            dataSpectrumNorm);
        writeCheckpoint(settings.CheckpointFile, checkpoint);
        lastCheckpointStep = iteration;
    end
    iteration = iteration + 1;
end

completedSteps = iteration - 1;
if strcmp(settings.Display, 'iter') && completedSteps > 0 && ...
        lastDisplayedStep ~= completedSteps
    displayIteration(completedSteps, difference, ...
        frequencyHistory(:, completedSteps + 1), lastDisplayedStep < 0);
    lastDisplayedStep = completedSteps;
end

[exitFlag, stopReason, stopMessage] = terminationStatus( ...
    completedSteps, difference, settings, ...
    outputStopRequested, timeLimitReached, toc(callTimer));

if ~isempty(settings.CheckpointFile) && ...
        lastCheckpointStep ~= completedSteps
    checkpoint = makeRestartState(settings, S, T, K, Alpha, weight, ...
        isRealInput, precision, nDC, initFreqType, initFreqMaximum, ...
        initialFrequencies, completedSteps, frequencyHistory, ...
        differenceHistory, coefficientSpectrum, phi, residual, ...
        dataSpectrumNorm);
    writeCheckpoint(settings.CheckpointFile, checkpoint);
    lastCheckpointStep = completedSteps;
end

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
info.InitialFrequencies = initialFrequencies;
info.DisplayInterval = settings.DisplayInterval;
info.TimeLimit = settings.TimeLimit;
info.ElapsedTime = toc(callTimer);
info.ExitFlag = exitFlag;
info.StopReason = stopReason;
info.Message = stopMessage;
info.CheckpointFile = settings.CheckpointFile;
if isempty(settings.CheckpointFile)
    info.LastCheckpointStep = [];
else
    info.LastCheckpointStep = lastCheckpointStep;
end
info.Iteration.steps = completedSteps;
info.Iteration.omega = frequencyHistory(:, 1:(completedSteps + 1));
info.Iteration.difference = differenceHistory(1:completedSteps);
info.Iteration.converged = exitFlag == 1;

restart = makeRestartState(settings, S, T, K, Alpha, weight, ...
    isRealInput, precision, nDC, initFreqType, initFreqMaximum, ...
    initialFrequencies, completedSteps, frequencyHistory, ...
    differenceHistory, coefficientSpectrum, phi, residual, ...
    dataSpectrumNorm);

if ~strcmp(settings.Display, 'off')
    fprintf('%s\n', stopMessage);
end
if ~isempty(settings.OutputFcn)
    progress = makeProgress(completedSteps, difference, ...
        frequencyHistory(:, completedSteps + 1), info.ElapsedTime);
    progress.stopReason = stopReason;
    invokeOutputFcn(settings.OutputFcn, progress, 'done');
end
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
addParameter(p, 'InitialFrequencies', []);
addParameter(p, 'Display', 'off');
addParameter(p, 'DisplayInterval', 20);
addParameter(p, 'OutputFcn', []);
addParameter(p, 'TimeLimit', inf);
addParameter(p, 'CheckpointFile', '');
addParameter(p, 'CheckpointInterval', 50);
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
addParameter(p, 'Device', state.Device);
addParameter(p, 'Display', state.Display);
addParameter(p, 'DisplayInterval', state.DisplayInterval);
addParameter(p, 'OutputFcn', []);
addParameter(p, 'TimeLimit', inf);
addParameter(p, 'CheckpointFile', '');
addParameter(p, 'CheckpointInterval', 50);
% Accepted for compatibility; restart state remains authoritative.
addParameter(p, 'InitFreqType', state.InitFreqType);
addParameter(p, 'InitFreqMaximum', state.InitFreqMaximum);
addParameter(p, 'FPPrecision', state.FPPrecision);
addParameter(p, 'nDC', state.nDC);
addParameter(p, 'InitialFrequencies', state.InitialFrequencies);
parse(p, varargin{:});

compatibility = p.Results;
compatibilityInitialFrequencies = validateInitialFrequencies( ...
    compatibility.InitialFrequencies, state.K, state.nDC);
if compatibility.InitFreqType ~= state.InitFreqType || ...
        compatibility.InitFreqMaximum ~= state.InitFreqMaximum || ...
        ~strcmpi(char(compatibility.FPPrecision), state.FPPrecision) || ...
        compatibility.nDC ~= state.nDC || ...
        ~isequal(cast(compatibilityInitialFrequencies, state.FPPrecision), ...
                 cast(state.InitialFrequencies, state.FPPrecision))
    error('rvmd:ImmutableRestartOption', ...
        ['Initialization frequencies, initialization settings, precision, ', ...
         'and nDC cannot change on restart.']);
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
settings.InitialFrequencies = state.InitialFrequencies;
settings.Display = validateChoice(p.Results.Display, ...
    {'off', 'final', 'iter'}, 'rvmd:InvalidDisplay', 'Display');
settings.DisplayInterval = validatePositiveInteger( ...
    p.Results.DisplayInterval, 'DisplayInterval');
settings.OutputFcn = validateOutputFcn(p.Results.OutputFcn);
settings.TimeLimit = validateTimeLimit(p.Results.TimeLimit);
settings.CheckpointFile = validateCheckpointFile(p.Results.CheckpointFile);
settings.CheckpointInterval = validatePositiveInteger( ...
    p.Results.CheckpointInterval, 'CheckpointInterval');
settings.Weight = state.weight;
end

function validateProblem(Q, K, Alpha)
if ~isnumeric(Q) || isempty(Q) || ~ismatrix(Q) || ...
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
settings.InitialFrequencies = validateInitialFrequencies( ...
    results.InitialFrequencies, K, settings.nDC);
settings.Display = validateChoice(results.Display, {'off', 'final', 'iter'}, ...
    'rvmd:InvalidDisplay', 'Display');
settings.DisplayInterval = validatePositiveInteger( ...
    results.DisplayInterval, 'DisplayInterval');
settings.OutputFcn = validateOutputFcn(results.OutputFcn);
settings.TimeLimit = validateTimeLimit(results.TimeLimit);
settings.CheckpointFile = validateCheckpointFile(results.CheckpointFile);
settings.CheckpointInterval = validatePositiveInteger( ...
    results.CheckpointInterval, 'CheckpointInterval');
end

function state = upgradeRestartState(input)
if ~isstruct(input) || ~isscalar(input)
    error('rvmd:InvalidRestart', ...
        'The restart input must be a scalar RVMD state struct.');
end
state = input;

if isfield(state, 'version')
    if ~isnumeric(state.version) || ~isscalar(state.version) || ...
            ~isfinite(state.version) || state.version < 1 || ...
            state.version ~= fix(state.version)
        error('rvmd:InvalidRestart', ...
            'The restart state has an invalid version value.');
    end
    sourceVersion = double(state.version);
else
    sourceVersion = 0;
end
if sourceVersion > currentRestartVersion()
    error('rvmd:UnsupportedRestartVersion', ...
        'Restart version %d is newer than this RVMD implementation.', ...
        sourceVersion);
end

required = {'S', 'T', 'K', 'alpha', 'weight', 'c_spec_n', 'phi_n'};
if sourceVersion < 4
    required{end + 1} = 'Q';
else
    required(end + (1:2)) = {'residual_n', 'dataSpectrumNorm'};
end
if ~all(isfield(state, required))
    error('rvmd:InvalidRestart', ...
        'The restart struct is missing required RVMD state fields.');
end

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
    if isa(state.phi_n, 'single') || isa(state.c_spec_n, 'single')
        state.FPPrecision = 'single';
    else
        state.FPPrecision = 'double';
    end
else
    state.FPPrecision = validateChoice(state.FPPrecision, ...
        {'single', 'double'}, 'rvmd:InvalidRestart', 'FPPrecision');
end
if ~isfield(state, 'isRealInput')
    if ~isfield(state, 'Q')
        error('rvmd:InvalidRestart', ...
            'The restart state does not identify whether its input was real.');
    end
    state.isRealInput = isreal(state.Q);
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
if ~isfield(state, 'DisplayInterval')
    state.DisplayInterval = 20;
end
if ~isfield(state, 'Device')
    state.Device = 'cpu';
end
if ~isfield(state, 'InitialFrequencies')
    state.InitialFrequencies = state.Iteration.omega(:, 1);
end

S = validatePositiveInteger(state.S, 'S');
T = validatePositiveInteger(state.T, 'T');
K = validatePositiveInteger(state.K, 'K');
if ~isnumeric(state.alpha) || ~isreal(state.alpha) || ...
        ~isscalar(state.alpha) || ~isfinite(state.alpha) || state.alpha < 0
    error('rvmd:InvalidRestart', ...
        'The saved bandwidth penalty must be finite and nonnegative.');
end
if ~islogical(state.isRealInput) || ~isscalar(state.isRealInput)
    error('rvmd:InvalidRestart', ...
        'The saved isRealInput flag must be a logical scalar.');
end
state.nDC = validateNonnegativeInteger(state.nDC, 'nDC');
if state.nDC > K
    error('rvmd:InvalidRestart', 'Saved nDC is greater than saved K.');
end
state.InitialFrequencies = validateInitialFrequencies( ...
    state.InitialFrequencies, K, state.nDC);

if ~isnumeric(state.Iteration.steps) || ...
        ~isscalar(state.Iteration.steps) || ...
        ~isfinite(state.Iteration.steps) || ...
        state.Iteration.steps < 0 || ...
        state.Iteration.steps ~= fix(state.Iteration.steps) || ...
        size(state.phi_n, 1) ~= S || size(state.phi_n, 2) ~= K || ...
        size(state.Iteration.omega, 1) ~= K || ...
        size(state.Iteration.omega, 2) < state.Iteration.steps + 1 || ...
        numel(state.Iteration.difference) < state.Iteration.steps
    error('rvmd:InvalidRestart', 'The restart struct is internally inconsistent.');
end
if isfield(state, 'Q') && ...
        (size(state.Q, 1) ~= S || size(state.Q, 2) ~= T)
    error('rvmd:InvalidRestart', ...
        'The saved Q dimensions do not match S and T.');
end
if state.isRealInput
    spectrumLength = T + 1;
else
    spectrumLength = 2 * T;
end
if ~isequal(size(state.c_spec_n), [spectrumLength, K]) || ...
        (isfield(state, 'residual_n') && ...
         ~isequal(size(state.residual_n), [S, spectrumLength]))
    error('rvmd:InvalidRestart', ...
        'The saved spectral state has inconsistent dimensions.');
end
if isfield(state, 'dataSpectrumNorm')
    state.dataSpectrumNorm = validateNonnegativeScalar( ...
        state.dataSpectrumNorm, 'DataSpectrumNorm');
end
state.Tolerance = validateNonnegativeScalar(state.Tolerance, 'Tolerance');
state.MaximumSteps = validatePositiveInteger( ...
    state.MaximumSteps, 'MaximumSteps');
if state.MaximumSteps < state.Iteration.steps
    error('rvmd:InvalidRestart', ...
        'Saved MaximumSteps is less than the completed iteration count.');
end
state.Display = validateChoice(state.Display, {'off', 'final', 'iter'}, ...
    'rvmd:InvalidRestart', 'Display');
state.DisplayInterval = validatePositiveInteger( ...
    state.DisplayInterval, 'DisplayInterval');
state.Device = validateChoice(state.Device, {'cpu', 'gpu'}, ...
    'rvmd:InvalidRestart', 'Device');
state.S = S;
state.T = T;
state.K = K;
end

function progress = makeProgress(step, difference, omega, elapsedTime)
progress.step = step;
progress.difference = double(difference);
progress.omega = gatherArray(omega);
progress.elapsedTime = double(elapsedTime);
progress.stopReason = '';
end

function stop = invokeOutputFcn(outputFcn, progress, phase)
stop = outputFcn(progress, phase);
if isempty(stop)
    stop = false;
end
if ~(islogical(stop) || isnumeric(stop)) || ~isscalar(stop) || ...
        ~isfinite(double(stop))
    error('rvmd:InvalidOutputFcnResult', ...
        'OutputFcn must return a finite logical or numeric scalar.');
end
stop = logical(stop);
end

function displayIteration(step, difference, omega, printHeader)
if printHeader
    fprintf(' Iteration    Difference      Center frequencies (cycles/sample)\n');
end
omega = sort(double(gatherArray(omega(:))).', 'ascend');
fprintf('%10d    %10.3e      %s\n', ...
    step, difference, mat2str(omega, 6));
end

function [exitFlag, stopReason, message] = terminationStatus( ...
        completedSteps, difference, settings, outputStop, timeStop, elapsedTime)
if completedSteps > 0 && difference <= settings.Tolerance
    exitFlag = 1;
    stopReason = 'converged';
    description = 'convergence tolerance satisfied';
elseif outputStop
    exitFlag = -1;
    stopReason = 'outputFunction';
    description = 'stopped by OutputFcn';
elseif timeStop
    exitFlag = -2;
    stopReason = 'timeLimit';
    description = 'TimeLimit reached';
else
    exitFlag = 0;
    stopReason = 'maximumSteps';
    description = 'MaximumSteps reached';
end
message = sprintf('RVMD %s after %d iterations (%.3g s).', ...
    description, completedSteps, elapsedTime);
end

function restart = makeRestartState(settings, S, T, K, Alpha, weight, ...
        isRealInput, precision, nDC, initFreqType, initFreqMaximum, ...
        initialFrequencies, completedSteps, frequencyHistory, ...
        differenceHistory, coefficientSpectrum, phi, residual, ...
        dataSpectrumNorm)
restart.version = currentRestartVersion();
restart.S = S;
restart.T = T;
restart.K = K;
restart.alpha = gatherArray(Alpha);
restart.weight = gatherArray(weight);
restart.Tolerance = settings.Tolerance;
restart.MaximumSteps = settings.MaximumSteps;
restart.InitFreqType = initFreqType;
restart.InitFreqMaximum = initFreqMaximum;
restart.InitialFrequencies = gatherArray(initialFrequencies);
restart.Device = settings.Device;
restart.FPPrecision = precision;
restart.nDC = nDC;
restart.isRealInput = isRealInput;
restart.Display = settings.Display;
restart.DisplayInterval = settings.DisplayInterval;
restart.Iteration.steps = completedSteps;
restart.Iteration.omega = gatherArray( ...
    frequencyHistory(:, 1:(completedSteps + 1)));
restart.Iteration.difference = gatherArray( ...
    differenceHistory(1:completedSteps));
if completedSteps == 0
    restart.Iteration.converged = false;
else
    restart.Iteration.converged = ...
        scalarValue(differenceHistory(completedSteps)) <= settings.Tolerance;
end
restart.c_spec_n = gatherArray(coefficientSpectrum);
restart.phi_n = gatherArray(phi);
restart.residual_n = gatherArray(residual);
restart.dataSpectrumNorm = double(dataSpectrumNorm);
end

function writeCheckpoint(checkpointFile, restart)
[checkpointDirectory, ~, ~] = fileparts(checkpointFile);
if ~isempty(checkpointDirectory) && exist(checkpointDirectory, 'dir') ~= 7
    error('rvmd:CheckpointWriteFailed', ...
        'Checkpoint directory does not exist: %s', checkpointDirectory);
end

temporaryFile = [checkpointFile, '.tmp'];
previousFile = [checkpointFile, '.prev'];
if exist(temporaryFile, 'file') == 2
    delete(temporaryFile);
end

try
    if exist('OCTAVE_VERSION', 'builtin') ~= 0
        save(temporaryFile, 'restart', '-v7');
    else
        save(temporaryFile, 'restart', '-v7.3');
    end
catch exception
    if exist(temporaryFile, 'file') == 2
        delete(temporaryFile);
    end
    error('rvmd:CheckpointWriteFailed', ...
        'Could not write checkpoint %s: %s', checkpointFile, exception.message);
end

if exist(previousFile, 'file') == 2
    delete(previousFile);
end
if exist(checkpointFile, 'file') == 2
    [moved, moveMessage] = movefile(checkpointFile, previousFile, 'f');
    if ~moved
        delete(temporaryFile);
        error('rvmd:CheckpointWriteFailed', ...
            'Could not preserve the previous checkpoint: %s', moveMessage);
    end
end
[moved, moveMessage] = movefile(temporaryFile, checkpointFile, 'f');
if ~moved
    if exist(previousFile, 'file') == 2 && ...
            exist(checkpointFile, 'file') ~= 2
        movefile(previousFile, checkpointFile, 'f');
    end
    error('rvmd:CheckpointWriteFailed', ...
        'Could not activate the new checkpoint: %s', moveMessage);
end
end

function version = currentRestartVersion()
version = 4;
end

function value = gatherArray(value)
if isa(value, 'gpuArray')
    value = gather(value);
end
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

function phi = normalizedFallback(phi, k, weight, isRealInput)
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

function frequencies = validateInitialFrequencies(frequencies, K, nDC)
if isempty(frequencies)
    frequencies = [];
    return
end
if ~isnumeric(frequencies) || ~isreal(frequencies) || ...
        ~isvector(frequencies) || numel(frequencies) ~= K || ...
        ~all(isfinite(frequencies(:))) || ...
        any(frequencies(:) < 0) || any(frequencies(:) > 0.5)
    error('rvmd:InvalidInitialFrequencies', ...
        'InitialFrequencies must be a finite K-vector in [0, 0.5].');
end
frequencies = double(frequencies(:));
if nDC > 0 && any(frequencies(1:nDC) ~= 0)
    error('rvmd:InvalidInitialFrequencies', ...
        'The first nDC initial frequencies must be zero.');
end
end

function outputFcn = validateOutputFcn(outputFcn)
if ~isempty(outputFcn) && ~isa(outputFcn, 'function_handle')
    error('rvmd:InvalidOutputFcn', ...
        'OutputFcn must be empty or a function handle.');
end
end

function value = validateTimeLimit(value)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        isnan(value) || value < 0
    error('rvmd:InvalidTimeLimit', ...
        'TimeLimit must be a nonnegative real scalar or Inf.');
end
value = double(value);
end

function checkpointFile = validateCheckpointFile(checkpointFile)
if ~isTextScalar(checkpointFile)
    error('rvmd:InvalidCheckpointFile', ...
        'CheckpointFile must be a text scalar.');
end
checkpointFile = char(checkpointFile);
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
