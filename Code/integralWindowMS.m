function [preStepMS, intStartMS, intDurMS] = integralWindowMS()
% return the analysis-wide limits on the kernel integration window
  preStepMS = 750;
  intStartMS = 0;
  intDurMS = 125;
end