function result = run_multiple_files()
myFolder = 'C:\Users\Owner\Documents\MATLAB\Synaptophysin_Image_Analysis';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.tif'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
A = zeros(2,56)
B = strings(1,56)
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  % Now do whatever you want with this file name,
  % such as reading it in as an image array with imread()
  X = imread(fullFileName);
  X = imresize(X,0.5);
  P = createMask_synap(X,fullFileName);% Display image.
  A(:,k) = P
  B(k) = fullFileName
  drawnow; % Force display to update immediately.
end
csvwrite('detected_pixels_synap.csv',A)
B
end