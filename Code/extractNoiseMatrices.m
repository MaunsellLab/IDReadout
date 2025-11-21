function [prefMat, probeMat, trialOutcomes] = extractNoiseMatrices(header, trials, stepTypes, sideType)
% extractNoiseMatrices  Compile selected coherence noise into a 2 x m x n matrix.
%
% We do partial bundling because the data are a bit complex.  The
% pref/probe trials for a given step direction (inc/dec) are always
% balanced because they come from the exact same trials.  But inc and dec
% trials might differ in number, so their matrices will have different
% sizes. Also, we might want to do analyses in the future that construct
% kernel based on left/right or other factors. So we've left this matrix
% formation fairly open, except for pairing the pref/prob matrices.
%
%   [prefMat, probMat] = extractNoiseMatrices(header, trials, changeType, incType)
%
%   INPUT:
%     header : cell array of header info from an IDKReadout data file
%     trials : cell array of trial info extracted from anIDKReadout data file
%     stepTypes : vector of valid stepTypes: coh decrement 1, increment 2 
%     sideType : side with step 0, RF step 1, or non-RF step 2
%
%   OUTPUT:
%     prefMat         : m x nTrials matrix for preferred noise
%     probeMat        : m x nTrials matrix for probe noise
%
%   Pref and probe noise are presented on the same trials, so their
%   matrices are always matched.  Inc/Dec or Change/NoChange might not be
%   matched
    
nTrials  = numel(trials);
if nTrials == 0
    error('extractNoiseMatrices:EmptyInput', 'Input "trials" is empty.');
end

  % ----- Find valid trials (trialEnd == 0 or 1, correct or wrong, for selected types) -----
  validIdx = [];
  trialOutcomes = [];  % storage for 0/1 correct/wrong trialEnd values
  for k = 1:nTrials
    tr = trials{k};
    if ~isfield(tr, 'trialEnd')
        continue;                       % skip if missing outcome
    end
    tCert = tr.trialCertify.data;       % trial certify (dropped frame monitor)
    tEnd = tr.trialEnd.data;            % end of trial code
    tSide = sideType == 0 || sideType == tr.trial.data.changeSide + 1;  % change index
    tStep = tr.trial.data.changeIndex + 1; % index for inc or dec step
    if tCert == 0 && ismember(tEnd, [0 1]) && tSide && ismember(tStep, stepTypes)
        validIdx(end + 1) = k;          %#ok<AGROW>
        trialOutcomes(end + 1) = tEnd;  %#ok<AGROW>
    end
  end
  nValid = numel(validIdx);
  if nValid == 0
      error('extractCohNoiseMatrices:NoValidTrials', 'No valid matching trials were found.');
  end

  % Preallocate output matrices
  msPerVFrame = 1000.0 / header.frameRateHz.data;
  m = round((header.preStepMS.data(1) + header.stepMS.data(1)) / msPerVFrame);  
  prefMat = nan(m, nValid);
  probeMat = nan(m, nValid);

  % ----- Process valid trials to extract coherence noise -----
  for k = 1:nValid
    tr = trials{validIdx(k)};
    % Take the noise from the change side if we're doing change side, or if
    % we're doing RF side and the change is there, or if we're doing the
    % opposite side and the change is there.
    if sideType == 0 || sideType == tr.trial.changeSide + 1   % process the side that changed
      prefCohsPC = tr.changePrefCohsPC.data(:);
      probeCohsPC = tr.changeProbeCohsPC.data(:);
      cohTimesMS = tr.changeTimesMS.data(:);
    else                                                      % process the other side
      prefCohsPC = tr.noChangePrefCohsPC.data(:);
      probeCohsPC = tr.noChangeProbeCohsPC.data(:);
      cohTimesMS = tr.noChangeTimesMS.data(:);
    end
  
    nTimes = length(cohTimesMS);
    for tIndex = 1:nTimes
      if tIndex == 1
        theVFrame = 1;
      else
        theVFrame = max(1, round(cohTimesMS(tIndex) / msPerVFrame));
      end
      if tIndex < nTimes
        nextVFrame = round(cohTimesMS(tIndex + 1) / msPerVFrame);
      else
        nextVFrame = m + 1;
      end
      for vFrame = theVFrame:nextVFrame - 1
          prefMat(vFrame, k) = prefCohsPC(tIndex);
          probeMat(vFrame, k) = probeCohsPC(tIndex);
      end
    end
  end
