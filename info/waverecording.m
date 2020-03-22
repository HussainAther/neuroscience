classdef wav_recording < recording
    methods 
        function obj = wav_recording(filename)
            obj = obj@recording(filename);
        end
        function r = sample_rate(obj)
            [data, r] = wavread(obj.filename);
        end
        function data = raw_data(obj)
            [data, r] = wavread(obj.filename);
        end
    end
end 

classdef pcm_recording < recording
    methods
        function obj = pcm_recording(filename, sample_rate)
            obj = obj@recording(filename);
            obj.sample_rate = sample_rate;
        end
        function r = sample_rate(obj)
            return obj.sample_rate;
        end
        function data = raw_data(obj)
            fid = fopen(obj.filename, "r");
            if fid == -1
                error("Unable to open file" + obj.filename);
                data = [];
                return
            else
                data = fread(fid, inf, "unit16 = > double", 0, "l");
                fclose(fid);
            end
        end
    end
end 

function dates = threshold_crossings(wav_filename, start_time)
    % Takes the WAV file wav_filename as input and the start time start_time
    % of the WAV recording.
    raw_data, sampling_raw = wavread(wav_filename);
    above = raw_data > threshold;
    threshold_crossings = diff(above) == 1;
    threshold_times = find(threshold_crossings);
    % threshold_times is the sample count since start_time
    % for each threshold crossing
