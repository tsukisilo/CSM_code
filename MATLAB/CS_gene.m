function [ Q ] = CS_gene( n )
%CS_GENE Summary of this function goes here
%   Detailed explanation goes here
Q = eye(n);
Q(:,2:n)=Q(:,1:n-1);
Q(n,1)=1;
Q(1,1)=0;
end

