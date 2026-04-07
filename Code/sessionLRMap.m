function lr = sessionLRMap(trials)
% sessionLRMap  Determine sessionwise mapping of RF/Opp onto Left/Right.
%
%   lr = sessionLRMap(trials)
%
% The mapping is fixed within a session. It is inferred from any valid trial.
%
% Convention:
%   Left if azimuthDeg < 0 OR (azimuthDeg == 0 AND elevationDeg >= 0)
%   Else Right
%
% changeSide coding:
%   0 = changing patch is RF
%   1 = changing patch is Opp
%
% Output:
%   lr.leftIsRF   true if Left patch is RF (and Right patch is Opp)
%   lr.rfAz       RF patch azimuth (for record/debugging)
%   lr.rfEl       RF patch elevation (for record/debugging)

  if nargin < 1 || isempty(trials)
    error('sessionLRMap:BadInput', 'trials is empty.');
  end

  t = [];
  for i = 1:numel(trials)
    if isempty(trials{i})
      continue;
    end
    ti = trials{i};
    if isfield(ti, 'changeDots') && isfield(ti, 'trial') && ...
       isfield(ti.changeDots, 'data') && isfield(ti.trial, 'data') && ...
       isfield(ti.changeDots.data, 'azimuthDeg') && ...
       isfield(ti.changeDots.data, 'elevationDeg') && ...
       isfield(ti.trial.data, 'changeSide')
      t = ti;
      break;
    end
  end

  if isempty(t)
    error('sessionLRMap:NoValidTrial', ...
      'Could not find a valid trial with changeDots and changeSide.');
  end

  az = t.changeDots.data.azimuthDeg;
  el = t.changeDots.data.elevationDeg;

  % changeSide: 0 = RF, 1 = Opp
  if t.trial.data.changeSide == 0
    rfAz = az;
    rfEl = el;
  else
    rfAz = -az;
    rfEl = -el;
  end

  leftIsRF = (rfAz < 0) || (rfAz == 0 && rfEl >= 0);

  lr = struct;
  lr.leftIsRF = leftIsRF;
  lr.rfAz = rfAz;
  lr.rfEl = rfEl;
end