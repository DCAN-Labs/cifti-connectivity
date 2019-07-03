function [tf, lthresh, uthresh, center] = isthisanoutlier(a,varargin)
% ISOUTLIER Find outliers in data
%   TF = ISOUTLIER(A) returns a logical array whose elements are true when 
%   an outlier is detected in the corresponding element. An outlier is an 
%   element that is greater than 3 scaled median absolute deviation (MAD) away 
%   from the median. The scaled MAD is defined as K*MEDIAN(ABS(A-MEDIAN(A))) 
%   where K is the scaling factor and is approximately 1.4826. If A is a 
%   matrix or a table, ISOUTLIER operates on each column separately. If A is
%   an N-D array, ISOUTLIER operates along the first array dimension
%   whose size does not equal 1.
%
%   TF = ISOUTLIER(A, METHOD) specifies the method to determine outliers, 
%   where METHOD is one of the following:
%     'median'    - Returns all elements more than 3 scaled MAD from the 
%                   median. This is the default.
%     'mean'      - Returns all elements more than 3 standard deviations 
%                   from the mean. This is also known as the three-sigma 
%                   rule, and is a fast but less robust method.
%     'quartiles' - Returns all elements more than 1.5 interquartile ranges 
%                   (IQR) above the upper quartile or below the lower quartile.
%                   This method makes few assumptions about data distribution,
%                   and is appropriate if A is not normally distributed.
%     'grubbs'    - Applies Grubbs' test for outliers, which is an iterative 
%                   method that removes one outlier per iteration until no
%                   more outliers are found. This method uses formal statistics
%                   of hypothesis testing and gives more objective reasoning backed
%                   by statistics behind its outlier identification. It assumes normal 
%                   distribution and may not be appropriate if A is not normal.
%     'gesd'      - Applies the generalized extreme Studentized deviate
%                   test for outliers. Like 'grubbs', this is another iterative
%                   method that removes one outlier per iteration. It may offer 
%                   improved performance over 'grubbs' when there are multiple
%                   outliers that mask one another.
%
%   TF = ISOUTLIER(A, MOVMETHOD, WL) uses a moving window method to determine 
%   contextual outliers instead of global outliers. Contextual outliers are 
%   outliers in the context of their neighborhood, and may not be the
%   maximum or minimum values in A. MOVMETHOD can be:
%     'movmedian' - Returns all elements more than 3 local scaled MAD from 
%                   the local median, over a sliding window of length WL.
%     'movmean'   - Returns all elements more than 3 local standard deviations 
%                   from the local mean, over a sliding window of length WL.
%   WL is the length of the moving window. It can either be a scalar or a
%   two-element vector, which specifies the number of elements before and
%   after the current element.
%
%   TF = ISOUTLIER(..., DIM) specifies the dimension to operate along.
%
%   TF = ISOUTLIER(..., 'ThresholdFactor', P) modifies the outlier detection 
%   thresholds by a factor P. For the 'grubbs' and 'gesd' methods, P is a 
%   scalar between 0 and 1. For all other methods, P is a nonnegative 
%   scalar. See the documentation for more information.
%
%   TF = ISOUTLIER(...,'SamplePoints',X) specifies the sample points 
%   X representing the location of the data in A, which is used by moving 
%   window methods. X must be a numeric or datetime vector, and must be 
%   sorted with unique elements. For example, X can specify time stamps for 
%   the data in A. By default, outliers uses data sampled uniformly at 
%   points X = [1 2 3 ... ].
%        
%   TF = ISOUTLIER(...,'DataVariables', DV) finds outliers only in the table 
%   variables specified by DV. The default is all table variables in A. DV 
%   must be a table variable name, a cell array of table variable names, a 
%   vector of table variable indices, a logical vector, or a function handle 
%   that returns a logical scalar (such as @isnumeric). TF has the same size as A.
%
%   TF = ISOUTLIER(..., 'MaxNumOutliers', MAXN) specifies the maximum number 
%   of outliers for the 'gesd' method only. The default is 10% of the number 
%   of elements. Set MAXN to a larger value to ensure it returns all outliers. 
%   Setting MAXN too large can reduce efficiency.
%
%   [TF, LTHRESH, UTHRESH, CENTER] = ISOUTLIER(...) also returns the
%   lower threshold, upper threshold, and the center value used by the 
%   outlier detection method.
%
%   Examples:
%      % Detect outliers in a data vector
%      x = [randn(1,50) 100 randn(1,49)];
%      tf = isoutlier(x);
%
%   Class support for input A:
%      float: double, single
%      table, timetable
%
%   See also FILLOUTLIERS, ISMISSING, FILLMISSING, RMMISSING, SMOOTHDATA

%   Copyright 2016 The MathWorks, Inc.

[method, wl, dim, p, sp, vars, maxoutliers] = parseinput(a, varargin);

xistable = matlab.internal.datatypes.istabular(a);

if xistable
    tf = false(size(a));
    if nargout > 1
        if ismember(method, {'movmedian', 'movmean'})
            % with moving methods, the thresholds and center have the same
            % size as input
            lthresh = a(:,vars);
        else
            % with other methods, thresholds and center has reduced
            % dimension along first dimension
            lthresh = a(1,vars);
        end
        uthresh = lthresh;
        center = lthresh;
    end
    for i = 1:length(vars)
        vari = a.(vars(i));
        if ~isfloat(vari)
            error(message('MATLAB:isoutlier:NonfloatTableVar',...
                a.Properties.VariableNames{vars(i)}, class(vari)));
        end
        [out, lt, ut, c] = locatealltheoutliers(vari, method, wl, p, ...
            sp, maxoutliers);
        tf(:,vars(i)) = any(out,2);
        if nargout > 1
            lthresh.(i) = lt;
            uthresh.(i) = ut;
            center.(i) = c;
        end
    end
else
    asparse = issparse(a);
    if dim > 1
        dims = 1:max(ndims(a),dim);
        dims(1) = dim;
        dims(dim) = 1;
        if asparse && dim > 2
            % permuting beyond second dimension not supported for sparse
            a = full(a);
        end
        a = permute(a, dims);
    end
    [tf, lthresh, uthresh, center] = locatealltheoutliers(a, method, ...
        wl, p, sp, maxoutliers);
    
    if dim > 1
        tf = ipermute(tf, dims);
        if asparse
            % explicitly convert to sparse. If dim > 2, we have converted
            % to full previously
            tf = sparse(tf);
        end
        if nargout > 1
            lthresh = ipermute(lthresh, dims);
            uthresh = ipermute(uthresh, dims);
            center = ipermute(center, dims);
            if asparse
                lthresh = sparse(lthresh);
                uthresh = sparse(uthresh);
                center = sparse(center);
            end
        end
    end
end

function [method, wl, dim, p, samplepoints, datavariables, maxoutliers] = ...
    parseinput(a, input)
method = 'median';
wl = [];
p = [];
dim = [];
samplepoints = [];
datavariables = [];
maxoutliers = [];
funcname = mfilename;

validateattributes(a,{'single','double','table','timetable'}, {'real'}, funcname, 'A', 1);
aistable = matlab.internal.datatypes.istabular(a);
if aistable
    datavariables = 1:width(a);
end

if ~isempty(input)
    i = 1;
    % parse methods and movmethod
    if ischar(input{i}) || isstring(input{i})
        str = validatestring(input{i},{'median', 'mean', 'quartiles', 'grubbs', ...
            'gesd', 'movmedian', 'movmean', 'SamplePoints', ...
            'DataVariables', 'ThresholdFactor', 'MaxNumOutliers'}, i+1);
        if ismember(str, {'median', 'mean', 'quartiles', 'grubbs','gesd'})
            % method
            method = str;
            i = i+1;
        elseif ismember(str, {'movmedian', 'movmean'})
            % movmethod
            method = str;
            if isscalar(input)
                error(message('MATLAB:isoutlier:MissingWindowLength',method));
            end
            wl = input{i+1};
            if (isnumeric(wl) && isreal(wl)) || islogical(wl) || isduration(wl) 
                if isscalar(wl)
                    if wl <= 0 || ~isfinite(wl) 
                        error(message('MATLAB:isoutlier:WindowLengthInvalidSizeOrClass'));
                    end
                elseif numel(wl) == 2
                    if any(wl < 0 | ~isfinite(wl)) 
                        error(message('MATLAB:isoutlier:WindowLengthInvalidSizeOrClass'));
                    end
                else
                    error(message('MATLAB:isoutlier:WindowLengthInvalidSizeOrClass'));
                end       
            else
                error(message('MATLAB:isoutlier:WindowLengthInvalidSizeOrClass'));
            end
            i = i+2;
        end
    end
    % parse dim
    if i <= length(input)
        if ~(ischar(input{i}) || isstring(input{i}))
            validateattributes(input{i},{'numeric'}, {'scalar', 'integer', 'positive'}, ...
                funcname, 'dim', i+1);
            dim = input{i};
            i = i+1;
        end
        
        % parse N-V pairs
        inputlen = length(input);
        if rem(inputlen - i + 1,2) ~= 0
            error(message('MATLAB:isoutlier:ArgNameValueMismatch'))
        end
        for i = i:2:inputlen
            name = validatestring(input{i}, {'SamplePoints', ...
                'DataVariables', 'ThresholdFactor', 'MaxNumOutliers'}, i+1);
            
            value = input{i+1};
            switch name
                case 'SamplePoints'
                    if istimetable(a)
                        error(message('MATLAB:isoutlier:SamplePointsTimeTable'));
                    end
                    if isfloat(value)
                        validateattributes(value,{'double','single'}, {'vector', 'increasing', 'finite', 'real'},...
                            funcname, 'SamplePoints', i+2)
                    elseif isdatetime(value) || isduration(value)
                        if ~(isvector(value) && issorted(value) &&  ...
                                length(unique(value))==length(value) && all(isfinite(value)))
                            error(message('MATLAB:isoutlier:InvalidSamplePoints'));
                        end
                    else
                        error(message('MATLAB:isoutlier:SamplePointsInvalidClass'));
                    end
                    samplepoints = value;
                case 'DataVariables'
                    if aistable
                        datavariables = unique(...
                            matlab.internal.math.checkDataVariables(a,value,'isoutlier'));
                    else
                        error(message('MATLAB:isoutlier:DataVariablesNonTable',class(a)));
                    end
                case 'ThresholdFactor'
                    validateattributes(value,{'numeric'}, {'real', 'scalar', ...
                        'nonnegative', 'nonnan'}, funcname, 'ThresholdFactor', i+2);
                    p = double(value);
                case 'MaxNumOutliers'
                    validateattributes(value,{'numeric'}, {'scalar', 'positive', ...
                        'integer'}, funcname, 'MaxNumOutliers', i+2);
                    maxoutliers = double(value);
            end
        end
    end
end
if isempty(p)  % default p
    switch method
        case {'median','mean','movmedian','movmean'}
            p = 3;
        case 'quartiles'
            p = 1.5;
        otherwise % grubbs, gesd
            p = 0.05;
    end
elseif ismember(method, {'grubbs', 'gesd'})
    if p > 1
        error(message('MATLAB:isoutlier:AlphaOutOfRange'));
    end
end

% dim
if isempty(dim)
    if aistable
        dim = 1;
    else
        dim = find(size(a) ~= 1,1);
        if isempty(dim)  % scalar x
            dim = 1;
        end
    end
elseif aistable && dim ~= 1
    error(message('MATLAB:isoutlier:TableDim'));
end

if ~isempty(maxoutliers)
    if ~strcmp(method, 'gesd')
        error(message('MATLAB:isoutlier:MaxNumOutliersGesdOnly'));
    elseif maxoutliers > size(a,dim)
        error(message('MATLAB:isoutlier:MaxNumOutliersTooLarge'));
    end
end

if ~isempty(samplepoints) && ~isequal(numel(samplepoints),size(a,dim))
    error(message('MATLAB:isoutlier:SamplePointsInvalidSize'));
end
if (isdatetime(samplepoints) || isduration(samplepoints)) && ...
        ~isempty(wl) && ~isduration(wl)
    error(message('MATLAB:isoutlier:SamplePointsNonDuration'));
end
if istimetable(a)
    samplepoints = a.Properties.RowTimes;
end

function [tf, lowerbound, upperbound, center, b] = locatealltheoutliers(a, method, ...
    wl, p, sp, maxoutliers, replace)
% LOCATEOUTLIERS Shared computation function for ISOUTLIER and FILLOUTLIERS

%   Copyright 2016 The MathWorks, Inc.
if islogical(method)
    % manual specification of outlier location
    tf = method;
    asiz = size(a);
    % propage sparsity of a using 'like'
    lowerbound = NaN([1 asiz(2:end)], 'like', a);
    upperbound = lowerbound;
    center = lowerbound;
else
    switch method
        case 'grubbs'
            asiz = size(a);
            ncols = prod(asiz(2:end));
            lowerbound = zeros([1 asiz(2:end)]);
            upperbound = lowerbound;
            center = lowerbound;
            aflat = a(:,:);
            
            tf = false(asiz);
            for i=1:ncols
                atemp = aflat(:,i);
                indvec = (i-1)*size(aflat,1)+1:i*size(aflat,1); % linear indices
                while true
                    n = length(atemp);
                    
                    center(i) = mean(atemp, 'omitnan');
                    astd = std(atemp, 'omitnan');
                    
                    adiff = abs(atemp - center(i));
                    [amax, loc] = max(adiff);
                    
                    t = datafuntinv(p/(2*n),n-2);
                    threshold = ((n-1)/sqrt(n))*abs(t)/sqrt(n-2+t^2);
                    
                    if amax/astd > threshold
                        atemp(loc) = [];
                        tf(indvec(loc)) = true;
                        indvec(loc) = [];
                    else
                        break;
                    end
                end
                lowerbound(i) = center(i) - astd*threshold;
                upperbound(i) = center(i) + astd*threshold;
            end
        case 'gesd'
            if isempty(maxoutliers)
                % Simply pick 10% of data size as maximum number of outliers
                maxoutliers = ceil(size(a,1)*0.1);
            end
            asiz = size(a);
            ncols = prod(asiz(2:end));
            lowerbound = NaN([1 asiz(2:end)]);
            upperbound = lowerbound;
            center = lowerbound;
            aflat = a(:,:);
            n = asiz(1);
            
            tf = false(asiz);
            if n > 0
                for j=1:ncols
                    indvec = (j-1)*size(aflat,1)+1:j*size(aflat,1); % linear indices
                    atemp = aflat(:,j);
                    amean = zeros(maxoutliers,1);
                    astd = zeros(maxoutliers,1);
                    lambda = zeros(maxoutliers,1);
                    R = zeros(maxoutliers,1);
                    Rloc = zeros(maxoutliers,1);
                    
                    for i = 1:maxoutliers
                        amean(i) = mean(atemp, 'omitnan');
                        astd(i) = std(atemp, 'omitnan');
                        [amax,loc] = max(abs(atemp - amean(i)));
                        R(i) = amax/astd(i);
                        atemp(loc) = [];
                        Rloc(i) = indvec(loc);
                        indvec(loc) = [];
                        
                        % compute lambda
                        pp = 1 - p / (2*(n-i+1));
                        t = datafuntinv(pp,n-i-1);
                        lambda(i) = (n-i)*t/sqrt((n-i-1+t.^2)*(n-i+1));
                    end
                    
                    lastindex = find(R > lambda, 1, 'last');
                    
                    tf(Rloc(1:lastindex)) = true;
                    if isempty(lastindex)
                        tindex = maxoutliers;
                    else
                        tindex = min(lastindex+1,maxoutliers);
                    end
                    center(j) = amean(tindex);
                    lowerbound(j) = amean(tindex) - astd(tindex)*lambda(tindex);
                    upperbound(j) = amean(tindex) + astd(tindex)*lambda(tindex);
                end
            end
            
        case 'median'
            madfactor = -1 /(sqrt(2)*erfcinv(3/2));  %~1.4826
            center = median(a,1,'omitnan');
            amad = madfactor*median(abs(a - center), 1, 'omitnan');
            
            lowerbound = center - p*amad;
            upperbound = center + p*amad;
        case 'mean'
            center = mean(a,1,'omitnan');
            astd = std(a,1,'omitnan');
            
            lowerbound = center - p*astd;
            upperbound = center + p*astd;
        case 'movmedian'
            madfactor = -1 /(sqrt(2)*erfcinv(3/2));  %~1.4826
            if isempty(sp)
                center = movmedian(a, wl, 1, 'omitnan');
                amovmad = madfactor*movmad(a, wl, 1, 'omitnan');
            else
                center = movmedian(a, wl, 1, 'omitnan', 'SamplePoints', sp);
                amovmad = madfactor*movmad(a, wl, 1, 'omitnan', 'SamplePoints', sp);
            end
            
            lowerbound = center - p*amovmad;
            upperbound = center + p*amovmad;
        case 'movmean'
            if isempty(sp)
                center = movmean(a, wl, 1, 'omitnan');
                amovstd = movstd(a, wl, 1, 'omitnan');
            else
                center = movmean(a, wl, 1, 'omitnan','SamplePoints',sp);
                amovstd = movstd(a, wl, 1, 'omitnan', 'SamplePoints', sp);
            end
            
            lowerbound = center - p*amovstd;
            upperbound = center + p*amovstd;
        otherwise %'quartiles'
            [aiqr, aquartiles] = datafuniqr(a);
            center = mean(aquartiles, 1, 'omitnan');  % used for replacement
            asiz = size(a);
            lquartile = reshape(aquartiles(1,:),[1,asiz(2:end)]);
            uquartile = reshape(aquartiles(2,:),[1,asiz(2:end)]);
            
            lowerbound = lquartile - p*aiqr;
            upperbound = uquartile + p*aiqr;
            
    end
    
    tf = (a > upperbound | a < lowerbound);
end

if nargout > 4
    % compute b
    b = a;
    if ischar(replace) || isstring(replace)
        switch replace
            case 'center'
                if ismember(method, {'movmedian', 'movmean'})
                    b(tf) = center(tf);
                else
                    b(tf) = center(ceil(find(tf)/size(a,1)));
                end
            case 'clip'
                b = min(max(b,lowerbound),upperbound);
            otherwise  % 'previous', 'next', 'nearest', 'linear', 'spline', 'pchip'
                % loop through columns
                b = b(:,:);  % flatten
                if isempty(sp)
                    sp = transpose(1:size(b,1));
                end
                isfloatsp = isfloat(sp);
                for colindex = 1:size(b,2)
                    % cast to full since GriddedInterpolant does not
                    % support sparse
                    bcol = full(b(:,colindex));                                                 
                    tfcol = tf(:,colindex);
                    numNonOutliers = sum(~tfcol);
                    if numNonOutliers > 1  % interpolation requires at least 2 data points
                        if isfloatsp
                            G = griddedInterpolant(sp(~tfcol),bcol(~tfcol),replace);
                            bcol(tfcol) = G(sp(tfcol)); % faster than interp1
                        else  % sp is datetime or duration
                            bcol(tfcol) = interp1(sp(~tfcol),bcol(~tfcol),sp(tfcol),replace,'extrap');
                        end
                    elseif numNonOutliers == 1  
                        % With one data point, we can replace for next,
                        % previous, or nearest. For the rest, do nothing.
                        nonOutlierIndex = find(~tfcol);
                        if strcmp(replace, 'nearest')
                            bcol(tfcol) = bcol(nonOutlierIndex);
                        elseif strcmp(replace, 'next')
                            bcol(1:nonOutlierIndex-1) = bcol(nonOutlierIndex);
                        elseif strcmp(replace, 'previous')
                            bcol(nonOutlierIndex+1:end) = bcol(nonOutlierIndex);
                        end
                    end
                    b(:,colindex) = bcol;
                end
                b = reshape(b, size(a));
        end
    else % isscalar(replace)
        b(tf) = replace;
    end
end