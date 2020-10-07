function [CRP, RESULTS, PARAMETERS]=mdcrqa(TS1,TS2,EMB,DEL,NORM,RAD,DLINE, VLINE, ZSCORE)
% mdcrqa
%  computes a cross-recurrence plot for two multi-dimensional time-series and
%  performs recurrence quantification:
%    [CRP, RESULTS, PARAMETERS]=mdcrqa_rev(TS1,TS2,EMB,DEL,NORM,RAD,DLINE, VLINE, ZSCORE)
%
%
% Inputs:
%
%  TS1 is a double-variable with each dimension of the to-be-correlated
%  signal as a row of numbers in a separate column.
%
%  TS1 is a double-variable with each dimension of the to-be-correlated
%  signal as a row of numbers in a separate column.
%
%  EMB is the number of embedding dimensions (i.e., EMB = 1 would be no
%  embedding via time-delayed surrogates, just using the provided number of
%  colums as dimensions.
%  The default value is EMB = 1.
%
%  DEL is the delay parameter used for time-delayed embedding (if EMB > 1).
%  The default value is DEL = 1.
%
%  NORM is the type of norm by with the phase-space is normalized. The
%  following norms are available:
%    'euc' - Euclidean distance norm
%    'max' - Maximum distance norm
%    'min' - Minimum distance norm
%    'non' - no normalization of phase-space
%  The default value is NORM = 'euc'.
%
%  DLINE is the minimum length of what is counted as a diagonal line (i.e.,
%  the minimum number diagonally adjacent data points) in order to compute
%  the diagonal line-based measures of the cross-recurrence plot.
%  The default value is DLINE = 2.
%
%  VLINE is the minimum length of what is counted as a vertical line (i.e.,
%  the minimum number vertically adjacent data points) in order to compute
%  the vertical line-based measures of the cross-recurrence plot.
%  The default value is DLINE = 2.
%
%  RAD is the threhold/radius size within points in phase-space are counted
%  as being recurrent.
%  The default value is RAD = 1.
%
%  ZSCORE indicats, whether the data (i.e., the different columns of TS1 and TS2,
%  being the different signals or dimensions of a signal) should be z-scored
%  before performing MdCRQA:
%    0 - no z-scoring
%    1 - z-score columns of TS1 and TS2
%  The default value is ZSCORE = 0.
%
%
% Outputs:
%
%  RP is a matrix holding the resulting recurrence plot.
%
%  RESULTS is a double-variable holding the following recurrence variables:
%    1.  Size of the CRP
%    2.  %REC  - percentage of recurrent points
%    3.  DNRLINES - number of diagonal lines
%    4.  %DET  - percentage of diagonally adjacent recurrent points
%    5.  MeanL - average length of adjacent recurrent points
%    6.  MaxL  - maximum length of diagonally adjacent recurrent points
%    7.  EntrL - Shannon entropy of distribution of diagonal lines
%    8.  VNRLINES - number of vertical lines
%    9.  %LAM  - percentage of vertically adjacent recurrent points
%    10.  MeanV - average length of diagonally adjacent recurrent points
%    11.  MaxV  - maximum length of vertically adjacent recurrent points
%    12. EntrV - Shannon entropy of distribution of vertical lines
%    
%  PARAMETERS is a cell-variable holding the employed parameter settings:
%    1. DIM
%    2. EMB
%    3. DEL
%    4. RAD
%    5. NORM
%    6. ZSCORE
%    
%
% Reference:
%
%  Wallot, S. (2018). Multidimensional Cross-Recurrence Quantification
%  Analysis (MdCRQA) - a method for quantifying correlation between
%  multivariate time-series. Multivariate Behavioral Research. Doi: 10.1080/00273171.2018.1512846
%
% Version:
%
% v1.0, 17. April 2018
% by Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics, Frankfurt, Germany
%
% The author does not give any warranty whatsoever on the reliability of
% the results obtained by this software.


if exist('TS1') % check whether input data has been specified - if not, throw error message
else
    error('No data has been specified for TS1.');
end

if exist('TS2') % check whether input data has been specified - if not, throw error message
else
    error('No data has been specified for TS2.');
end

testData = size(TS1)-size(TS2);

if testData(1) == 0 % check wehter input signals match in dimensionality and data points
else
    error('Input data (TS1 and TS2) do not match in terms of dimensionality (number of columns).');
end

if testData(2) == 0 % check wehter input signals match in dimensionality and data points
else
    error('Input data (TS1 and TS2) do not match in terms of number of data points (number of rows).');
end

clear testData

if exist('EMB') % check whether EMB has been specified - if not, set EMB = 1 (no surrogate embedding)
else
    EMB=1;
end

if exist('DEL') % check whether DEL has been specified - if not, set DEL = 1 (no delaying)
else
    DEL=1;
end

if length(TS1) > (EMB-1)*DEL
else
    error('Insufficient number of data points to embedd time-series (number of data points < (EMB-1)*DEL).');
end

if exist('NORM')  % check whether NORM has been specified - if not, set NORM = 'euc'; if yes, check whether specification is appropriate
    if NORM == 'euc' | NORM == 'min' | NORM == 'max' | NORM == 'non'
    else
        error('No appropriate norm parameter specified.');
    end
else
    NORM == 'euc';
end

if exist('RAD') % check whether RAD has been specified - if not, set arbitarily RAD = 1
else
    RAD=1;
end

if exist('ZSCORE') % check whether ZSCORE has been specified - if not, don't zscore
    if ZSCORE == 0
    else
    TS1 = zscore(TS1);
    TS2 = zscore(TS2);
    end
else
end

DIM = size(TS1); % compute dimensionality of input signals
DIM = DIM(2);

if EMB > 1 % if EMB > 1, perform time-delayed embbedding
    for i = 1:EMB
        tempTS1(1:length(TS1)-(EMB-1)*DEL,1+DIM*(i-1):DIM*i)=TS1(1+(i-1)*DEL:length(TS1)-(EMB-i)*DEL,:);
    end
    TS1=tempTS1;
    clear tempTS1
        for i = 1:EMB
        tempTS2(1:length(TS2)-(EMB-1)*DEL,1+DIM*(i-1):DIM*i)=TS2(1+(i-1)*DEL:length(TS2)-(EMB-i)*DEL,:);
    end
    TS2=tempTS2;
    clear tempTS2
else
end

PARAMETERS={DIM,EMB,DEL,NORM,RAD,DLINE,VLINE,ZSCORE}; % store parameters

a=pdist2(TS1,TS2); % create distance matrix and cross-recurrence matirx
a=abs(a)*-1;
if NORM == 'euc'
    b=mean(a(a<0));
    a=a/abs(b);
elseif NORM == 'min'
    b=max(a(a<0));
    a=a/abs(b);
elseif NORM == 'max'
    b=min(a(a<0));
    a=a/abs(b);
elseif NORM == 'non'
else
end
a=a+RAD;
a(a >= 0) = 1;
a(a < 0) = 0;
diag_hist=[];
vertical_hist=[];
RESULTS(1,1)=length(a); % calculate size of recurrence plot
RESULTS(1,2)=100*(sum(sum(a)))/(length(TS1)^2); % calculate percent recurrence
if RESULTS(1,2) > 0 % check whether percent recurrence is bigger 0
for i = -(length(TS1)-1):length(TS1)-1 % caluculate diagonal line distribution
    c=diag(a,i);
    d=bwlabel(c,8);
    d=tabulate(d);
    if d(1,1)==0
        d=d(2:end,2);
    else
        d=d(2);
    end
    diag_hist(length(diag_hist)+1:length(diag_hist)+length(d))=d;
end
diag_hist=diag_hist(diag_hist>=DLINE);
if isempty(diag_hist)
    diag_hist=0;
else
end

for i=1:length(TS1) % calculate vertical line distribution
    c=(a(:,i));
    v=bwlabel(c,8);
    v=tabulate(v);
    if v(1,1)==0
        v=v(2:end,2);
    else
        v=v(2);
    end
    vertical_hist(length(vertical_hist)+1:length(vertical_hist)+length(v))=v;
end
vertical_hist=vertical_hist(vertical_hist>=VLINE);
if isempty(vertical_hist)
    vertical_hist=0;
else
end

RESULTS(1,3)=length(diag_hist);
RESULTS(1,4)=100*sum(diag_hist(diag_hist>1))/sum(diag_hist); % calculate recurrence variables
RESULTS(1,5)=mean(diag_hist(diag_hist>1));
RESULTS(1,6)=max(diag_hist);
[count,bin]=hist(diag_hist(diag_hist>1),min(diag_hist(diag_hist>1)):max(diag_hist));
total=sum(count);
p=count./total;
del=find(count==0); p(del)=[];
RESULTS(1,7)=-1*sum(p.*log2(p));
RESULTS(1,8)=length(vertical_hist);
RESULTS(1,9)=100*sum(vertical_hist(vertical_hist>1))/sum(vertical_hist);
RESULTS(1,10)=mean(vertical_hist(vertical_hist>1));
RESULTS(1,11)=max(vertical_hist);
[count,bin]=hist(vertical_hist(vertical_hist>1),min(vertical_hist(vertical_hist>1)):max(vertical_hist));
total=sum(count);
p=count./total;
del=find(count==0); p(del)=[];
RESULTS(1,12)=-1*sum(p.*log2(p));
else
    RESULTS(1,3:12)=NaN;
end

CRP=imrotate(a,90); % format recurrence plot
end