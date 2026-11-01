function value = poolSplitSignalResponse(response,directionsDeg,modelName)
% POOLSPLITSIGNALRESPONSE  Competing readouts of the same noisy MT landscape.

switch lower(modelName)
    case 'max'
        value=max(response,[],2);
    case 'vector'
        basis=single(exp(1i*deg2rad(directionsDeg(:))));
        value=abs(response*basis)/single(numel(directionsDeg));
    case 'energy'
        value=mean(response.^2,2);
    otherwise
        error('Unknown pooling model: %s',modelName);
end
end
