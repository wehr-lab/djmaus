function Pupil_Dilation(dir, video_file)
% loads mp4 file
% requires user input to find pupil
% ira 8.31.17
commandStr = sprintf('python E:\\djmaus-data\\pupil_dilation.py "%s" "%s"', dir, video_file);

[status, commandOut] =system(commandStr);
if status ==0
    fprintf('loaded video')
elseif status == 2
    fprintf('failed to load video, could not find it')
end
