function saveMovieCube(movieCube, filename)
%  saveMovieCube(movieCube, filename)
%  This function saves the movie cube as a raw YUV-format 
%  video file for later conversion. One second of black is
%  added at the beginning to avoid start-up issues
%
%
fprintf('Saving your movie cube as %s\n', filename);

fprintf('This will take a while, please be patient\n');


fid = fopen(filename, 'w');
if fid < 0
  fprintf('Error opening file!\n');
  return;
end;

 
s = size(movieCube);


fprintf('Writing movie frames to raw file...\n');
uv = ones(s(1)/2, s(2)/2, 'uint8'); 
uv(:, :) = 128; 

blackf = zeros(s(1), s(2), 'uint8'); 
for frame = 1:30
  fwrite(fid, blackf', 'uchar'); 
     fwrite(fid, uv, 'uchar'); 
     fwrite(fid, uv, 'uchar'); 
end;   

for frame = 1:s(3) 
     nframe = movieCube(:, :, frame);
     fwrite(fid, nframe', 'uchar');
     % now write the constant luminance data
     fwrite(fid, uv, 'uchar'); 
     fwrite(fid, uv, 'uchar'); 
     
end


fclose(fid);
