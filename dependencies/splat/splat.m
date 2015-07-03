function [ y ] = splat( x )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
	y = mat2cell(x, size(x, 1), ones(1, size(x,2)));
end

