
classdef metrics
    
    methods (Static)
    
    function [ d ] = distancy(u_train, v_train, tau)
        % returns the distance between two points
        % tau is a timescale expressing the precision of spike times
        
        d = square_term(u_train, tau) + square_term(v_train, tau) - 2*cross_term(u_train, v_train, tau);
    end
   
    
    function [ square ] = square_term(u_train, tau)
        % returns the squre term of a spike train
        % tau is a timescale expressing the precision of spike times
        
        square = 0;
        u_len = length(u_train);
        
        for ui = 1:u_len
            for uj = 1:u_len
                square = square + exp(-abs(u_train(ui) - u_train(uj)) / tau);
            end
        end
    end
    

    function [ cross ] = cross_term(u_train, v_train, tau)
        % returns the crosss term between spike train u and v
        % tau is a timescale expressing the precision of spike times
        
        cross = 0;
        u_len = length(u_train);
        v_len = length(v_train);
        
        for ui = 1:u_len
            for vj = 1:v_len
                cross = cross + exp(-abs(u_train(ui) - v_train(vj)) / tau);
            end
        end
    end
    end
end


