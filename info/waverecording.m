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
