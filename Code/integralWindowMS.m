function [preStepMS, intStartMS, intDurMS] = integralWindowMS()
% return the analysis-wide limits on the kernel integration window
  preStepMS = 750;
  intStartMS = 50;
  intDurMS = 50;
end