function allocs

steps = [
  3, 8, 13, 20, 35;
  3, 8, 11, 16, 30
  3, 7, 12, 18, 30;
  3, 6, 9, 13, 20,
  ];

allocs = [0.500, 0.625, 0.750, 1.000];

nSteps = length(steps);
nAllocs = length(allocs);

showAllocs(steps, allocs, 0.0);
showAllocs(steps, allocs, 7.071);

end

function showAllocs(steps, allocs, noise)

fprintf('Noise: %.1f%%:\n', noise);
for aIndex = 1:length(allocs)
  alloc = allocs(aIndex);
  fprintf('Alloc %.3f/%.3f: ', alloc, 1.0 - alloc);
  for sIndex = 1:length(steps)
    allocStep = steps(aIndex, sIndex) * alloc + noise;
    otherStep = steps(aIndex, sIndex) * (1.0 - alloc) + noise;
    fprintf('%8.3f', allocStep / (allocStep + otherStep));
  end
  fprintf('\n');
end
fprintf('\n');
end