function figureHandle = rvmdplot(info, varargin)
%RVMDPLOT Plot RVMD convergence diagnostics.
%   RVMDPLOT(INFO) creates a figure containing INFO.Iteration.difference and
%   the unsorted INFO.Iteration.omega trajectories returned by RVMD.
%
%   RVMDPLOT(INFO, 'SampleRate', FS) scales frequencies from cycles per
%   sample by a finite positive scalar FS. Default: 1. When FS is 1, the
%   frequency-axis label is cycles/sample.
%
%   FIGUREHANDLE = RVMDPLOT(...) returns the created figure handle.

p = inputParser;
p.FunctionName = 'rvmdplot';
addRequired(p, 'info');
addParameter(p, 'SampleRate', 1);
parse(p, info, varargin{:});

validateInfo(info);
sampleRate = p.Results.SampleRate;
if ~isnumeric(sampleRate) || ~isreal(sampleRate) || ...
        ~isscalar(sampleRate) || ~isfinite(sampleRate) || sampleRate <= 0
    error('rvmdplot:InvalidSampleRate', ...
        'SampleRate must be a finite positive real scalar.');
end

steps = info.Iteration.steps;
difference = double(info.Iteration.difference(:));
omega = double(info.Iteration.omega) * double(sampleRate);

figureHandle = figure('Color', 'white');

subplot(2, 1, 1);
if isempty(difference)
    plot(nan, nan);
else
    semilogy(1:steps, max(difference, realmin), 'LineWidth', 1.2);
end
grid on;
xlabel('Iteration');
ylabel('Difference');
title(sprintf('RVMD convergence: %s', info.StopReason), ...
    'Interpreter', 'none');

subplot(2, 1, 2);
plot(0:steps, omega.', 'LineWidth', 1.1);
grid on;
xlabel('Iteration');
if sampleRate == 1
    ylabel('Center frequency (cycles/sample)');
else
    ylabel('Center frequency');
end
end

function validateInfo(info)
required = {'Iteration', 'StopReason'};
if ~isstruct(info) || ~all(isfield(info, required)) || ...
        ~isstruct(info.Iteration) || ...
        ~all(isfield(info.Iteration, {'steps', 'difference', 'omega'}))
    error('rvmdplot:InvalidInfo', ...
        'Input must be an info struct returned by RVMD.');
end
end
