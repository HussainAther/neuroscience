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
