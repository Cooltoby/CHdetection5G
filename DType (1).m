% DType.m
%========

classdef DType < int32
    enumeration
        invalid (-3)        
        mean(-2)
        deterministic(-1)
        numeric(0)
        exp(1)
        exponential(1)
        gauss(2)
        gaussian(2)
        gamma(3)
        uniform(4)
        beta(5)
        fishertippet(6)
        gev(6)
        waitingDet(7)
    end
end