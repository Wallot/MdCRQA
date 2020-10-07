function [lag, REC]=dcrp(CRP,lags,ang)
% dcrp
%  computes the diagonal cross-recurrence profile for a given
%  cross-recurrence plot:
%    [lag,REC]=dcrp(CRP,lags,ang)
%
%
% Inputs:
%
%  CRP is a binary (cross-)recurrence matrix
%
%  lags is the number of lags around the diagonal of that matrix for which
%  percent recurrence should be computed
%  
%  ang is the angle with which the cross-recurrence matrix is rotated.
%  If the input matrix is rotated so that the main diagonal (cross-recurrence at lag0) runs
%  from lower-left to upper-right, the matrix needs to be rotated by 90 degrees (i.e., ang = 1).
%  Otherwise, the input matrix does not need to be rotated (i.e., ang = 0).
%  The plots from the mdcrqa.m-function are rotated by 90 degrees. Hence, ang = 1.
%
%
% Outputs:
%
%  REC = percent recurrence at a respective diagonal around the central
%  diagonal of the recurrence matrix
%
%  lag = lag index for REC
%
%
% Reference:
%
%  Wallot, S. (2017). Multidimensional Cross-Recurrence Quantification
%  Analysis (MdCRQA) - a method for quantifying correlation between
%  multivariate time-series. ???
%
%
% Version:
%
% v1.0, 04. October 2017
% by Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics, Frankfurt, Germany

if exist('CRP') % check whether input data has been specified - if not, throw error message
else
    error('Recurrence matrix is missing.');
end

if exist('lags') % check whether lags has been specified - if not, throw error message
else
    error('lags has not been specified.');
end

if exist('ang') % check whether ang has been specified - if not, set to default value
    if ang == 0
    else
    ang = 1;
    end
end

if sum(lags < size(CRP)) == 2 % check whether sufficient data is availalbe
else
    error('lag is bigger than number of diagonal lines on the plot.')
end

if ang == 1
CRP = imrotate(CRP, 90);
else
end

j = 1;
for i = -lags:1:lags % caluculate diagonal recurrences at lag i
    REC(j)=100*sum(diag(CRP,i))/length(diag(CRP,i));
    j = j+1;
end

lag = -lags:1:lags; % compute lag-vector

end