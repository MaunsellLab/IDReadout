function [preStepMS, intStartMS, intDurMS] = integralWindowMS()
% return the analysis-wide limits on the kernel integration window
  preStepMS = 750;
  intStartMS = 25;
  intDurMS = 100;
end