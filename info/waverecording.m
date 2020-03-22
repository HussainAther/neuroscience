classdef wav_recording < recording
    methods 
        function obj = wav_recording(filename)
            obj = obj@recording(filename);
        end
      
