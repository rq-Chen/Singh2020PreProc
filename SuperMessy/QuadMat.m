function [Out] = QuadMat(X,P)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% computes x^T P x for multiple X's
PX=P*X;

Out=sum(X.*PX,1);
%for i=1:size(X,2)
%    Out(i)=X(:,i)'*PX(:,i);
%end
    
    
end

