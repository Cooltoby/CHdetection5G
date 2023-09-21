
% Part of the Distribution Class for Queuing Theory 
%==================================================

classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) Distribution
    
    
    % -- add support parameter
    %
    % ------------------------------------------------------------------- %
    % PROPERTIES -------------------------------------------------------- %
    % ------------------------------------------------------------------- %    
    properties %(SetAccess = private, GetAccess = public)
        Type    % DType
        Params  % Vector of doubles
        ParamNames % Cell Vector containing the names of the parameters
%         Support % support of the pdf
        SupportLeft % support of the pdf
        SupportRight % support of the pdf
        DetermX = []
        DetermY = []
        PdfX    %
        PdfY    %
        CdfX    %
        CdfY    %
%         MeanVal % 
%         VarVal  %
        PdfTol; 
        PdfGridNumber; 
        isdefault = true;
        NonNormalized = false;
        NumberOfSamplesForCreation = [nan; nan]; % if this was created from samples.
        ID
        created
    end
    % ------------------------------------------------------------------- %
    % CONSTANT PROPERTIES ----------------------------------------------- %
    % ------------------------------------------------------------------- %
    properties (Constant)
        DiscTol1 = 1e-3; % tolerance to form clusters
        DiscTol2 = 1e-9; % tolerance to remove single discrete values
        DistributionWarningLvl = 1; % 0.. no warnings. 1.. important warnings. 2.. all warnings.
        DistributionWarningTol = 1e-1;
        EnableDistributionWarnings = true;        
        DebugMode = false;
        IntegrationPoints = 2;
        DefaultPdfTol = 1e-9; 
%         DefaultPdfGridNumber = 1e3; 
%         DefaultPdfX = linspace(-1000,1000,20001);        
%         DefaultPdfX = linspace(-100,100,20001);        
%         DefaultPdfX = linspace(0,200,20001);                
%         DefaultPdfX = linspace(0,200,200001);     
%         DefaultPdfX = linspace(0,400,40001);      
        
%         DomainChanged = true;
% %         DomainChanged = false;
% %         DefaultPdfX = linspace(0,200,40001); 
% %         DefaultPdfX = linspace(-200,200,40001); 
% %         DefaultPdfX = linspace(-200,200,80001); 
% 
% %         DefaultPdfX = linspace(-200,200,160001); 
% 
% %         DefaultPdfX = linspace(-200,600,160001); 
% %         DefaultPdfX = linspace(-1000,2000,210001); 
%         DefaultPdfX = linspace(-1000,2000,60001); % THIS WAS THE LAST DEFAULT;
        
%         DefaultPdfX = linspace(-400,1000,160001); 
        
%         DefaultPdfX = linspace(-200,200,16001); 
        
%         DefaultPdfX = linspace(0,1000,40001);     
%         DefaultPdfX = linspace(0,4000,40001);      
% % % % %         DefaultPdfX = linspace(0,8000,40001);      
        
%         DefaultPdfX = linspace(0,4000,4001);      
%         DefaultPdfX = -1000:0.1:1000;        
    end
   
    % ------------------------------------------------------------------- %
    % ORDINARY METHODS -------------------------------------------------- %
    % ------------------------------------------------------------------- %
    methods
        %%
        function obj = Distribution(type, varargin)
            % Constructor for a Distribution object
            % type is one of the Types in DType
            % depending on type, different parameters can be given as input
            
            obj.ID = Distribution.currentID();
            obj.created = now;
            
            if Distribution.DebugMode
                fprintf('DEBUG: Distribution(%d)\n', obj.ID);
            end
            
            if nargin == 0
%                 warning('Invalid Distribution created');
                obj.Type = DType.invalid;
                return
            end
            
%             if Distribution.EnableDistributionWarnings && Distribution.DomainChanged
            if Distribution.DistributionWarningLvl > 1 && Distribution.DomainChanged
                warning('Different Default Domain!');
            end
                                    
            obj.PdfTol = Distribution.DefaultPdfTol; 
            X = Distribution.DefaultPdfX;
%             obj.PdfGridNumber = Distribution.DefaultPdfGridNumber;
%             obj.Support = [-Inf, +Inf];
            obj.SupportLeft = -Inf;
            obj.SupportRight= +Inf;
            obj.Type = DType(type);
            switch (obj.Type)
                case DType.exp
                    assert(nargin == 2);
                    lambda = varargin{1};
                    % assert(lambda > 0); % allow "negative exponentials"
                    assert(lambda ~= 0);
                    obj.Params = lambda;
                    obj.ParamNames = {'lambda (mean arrival rate)'};
                    if lambda > 0
                        obj.SupportLeft = 0;
                        obj.SupportRight= +Inf;
                    else
                        obj.SupportLeft = -Inf;
                        obj.SupportRight= 0;
                    end
                case DType.gamma
                    assert(nargin == 3);
                    p = varargin{1};
                    b = varargin{2};
                    assert(p > 0 && b > 0)
                    obj.Params = [p b];
                    obj.ParamNames = {'p (shape)', 'b (rate)'};
                    obj.SupportLeft = 0;
                    obj.SupportRight= +Inf;
                case DType.beta
                    assert(nargin >= 3);
                    alp = varargin{1};
                    bet = varargin{2};                    
                    assert(alp > 0 && bet > 0)
                    if nargin >= 4
                        a = varargin{3};
                        b = varargin{4};
                    else
                        a = 0;
                        b = 1;
                    end
                    assert(b>a);
                    obj.Params = [alp bet a b];
                    obj.ParamNames = {'alpha (shape)', 'beta (shape)', 'a (lower bound)', 'b (upper bound)'};
                    obj.SupportLeft = a;
                    obj.SupportRight= b;    
                case DType.gauss
                    m = 0; v = 1;
                    if nargin > 1
                        m = varargin{1};
                    end
                    if nargin > 2
                        v = varargin{2};
                    end
                    assert(v > 0);
                    obj.Params = [m, v];
                    obj.ParamNames = {'mu (mean)', 'sigma^2 (variance)'};
%                     obj.Support = [-Inf, + Inf];
%                     obj.SupportLeft = -Inf;
%                     obj.SupportRight= +Inf;
                case DType.mean
                    assert(nargin == 2);
                    obj.Params = varargin{1};
                    obj.ParamNames = {'mean'};
                case DType.deterministic
                    assert(nargin >= 2);
                    if nargin == 2
                        obj.DetermX = varargin{1};
                        obj.DetermY = 1;
                    else
                        obj.DetermX = varargin{1};
                        obj.DetermY = varargin{2};
                        if nargin > 3
                            obj.NonNormalized = varargin{3};
                        end
                        if ~obj.NonNormalized && Distribution.EnableDistributionWarnings && abs(sum(obj.DetermY) - 1) > Distribution.DistributionWarningTol
                            warning('Probabilities do not sum up to 1. (err=%g)', abs(sum(obj.DetermY) - 1));
                        end
                    end
                    obj.Params = [];
                    % compatible support with the domain:
                    spl =  X( find( X <= min([obj.DetermX]), 1, 'first') );
                    spr =  X( find( X >= max([obj.DetermX]), 1, 'last') );
                    if isempty(spl)
                        spl = min([obj.DetermX]);
                    end
                    if isempty(spr)
                        spr = max([obj.DetermX]);
                    end
                    obj.SupportLeft = spl;
                    obj.SupportRight= spr;
                    %
                case DType.numeric
                    assert(nargin >= 3);
                    obj.Params = [];
                    obj.PdfX = varargin{1};
                    if isempty(obj.PdfX)
                        obj.PdfX = Distribution.DefaultPdfX;
                    end
                    obj.PdfY = varargin{2};
                    assert(length(obj.PdfX) == length(obj.PdfY));
%                     obj.Support = [obj.PdfX(1), obj.PdfX(end)];
                    if nargin > 4
                        %deterministic parts:
                        obj.DetermX = varargin{3};
                        obj.DetermY = varargin{4};
                    end
                    if ~isempty(obj.DetermX)
                        X = obj.PdfX;
                        % compatible support with the domain:
                        spl =  X( find( X <= min([obj.DetermX]), 1, 'first') );
                        spr =  X( find( X >= max([obj.DetermX]), 1, 'last') );
                        if isempty(spl)
                            spl = min([obj.DetermX]);
                        end
                        if isempty(spr)
                            spr = max([obj.DetermX]);
                        end
                        obj.SupportLeft = spl;
                        obj.SupportRight= spr;
                        %
                    else
                        spl = []; spr = [];
                    end
                    if nargin > 5
                        obj.NonNormalized = varargin{5};
                    end
                    obj.SupportLeft = min([obj.PdfX(1) spl]);
                    obj.SupportRight= max([obj.PdfX(end) spr]);
                    obj.isdefault = isequal(obj.PdfX, Distribution.DefaultPdfX);
                case DType.uniform                    
                    assert(nargin == 3);
                    l = varargin{1};
                    r = varargin{2};
                    obj.Params = [l r];
                    obj.ParamNames = {'l (lower bound)', 'r (upper bound)'};
                    obj.SupportLeft = l;
                    obj.SupportRight= r;
                case DType.fishertippet
                    assert(nargin >= 3);
                    mu = varargin{1};
                    sigma=varargin{2};
                    if nargin > 3
                        xi = varargin{3};
%                         warning('to be implemented!');
                    else
                        xi = 0;
                    end
                    obj.Params = [mu, sigma, xi];                    
                    obj.ParamNames = {'mu (location)', 'sigma (scale)', 'xi (shape)'};
                    if xi > 0 
                        obj.SupportLeft = mu - sigma / xi;
                        obj.SupportRight= +inf;                        
                    elseif xi < 0 
                        obj.SupportLeft = -inf;
                        obj.SupportRight= mu - sigma / xi;                        
                    else
                        obj.SupportLeft = -inf;
                        obj.SupportRight= +inf;
                    end
                case DType.waitingDet
                    assert(nargin >= 3);
                    lambda = varargin{1};
                    servicetime =varargin{2};
                    if nargin > 3
                        c = varargin{3};
%                         warning('to be implemented!');
                    else
                        c = 1;
                    end                    
                    obj.Params = [lambda servicetime c];                    
                    obj.ParamNames = {'lambda (mean arrival rate)', 'D (service time)', 'c (number of servers)'};
                    obj.SupportLeft = 0;
                    obj.SupportRight= Inf;
                case DType.invalid
%                     warning('invalid distribution!')
                    return
                otherwise
                    error('Unknown Distribution Type!');
            end
            % unify deterministic parts:
            [C, ia, ic] = unique(obj.DetermX);
            obj.DetermX = C;
            determY = zeros(size(C));
            for i=1:length(ic)
                determY(ic(i)) = determY(ic(i)) + obj.DetermY(i);
            end
            obj.DetermY = determY;
            % remove zeros from deterministic values:
            isz = obj.DetermY == 0;
            obj.DetermX(isz) = [];
            obj.DetermY(isz) = [];
            %
            obj = obj.updatePdf();
            obj = obj.updateCdf();
        end        
        %% --------------------------------------------------------------- %
        function obj = updatePdf(obj, X)
            if Distribution.DebugMode
                fprintf('DEBUG: updatePdf([%s])\n', sprintf('%d ',  obj.ID));
            end
            
             % Handling of arrays of Distributions:
            if numel(obj) > 1
                res = Distribution.initDist(size(obj));
                for oi=1:numel(res)
                    if nargin < 2
                        res(oi) = obj(oi).updatePdf();
                    else
                        res(oi) = obj(oi).updatePdf(X);
                    end
                end
                obj = res;
                return;
            end
            %
            
            if obj.Type == DType.invalid 
%                 obj = Distribution(); %return invalid
                return;
            end
            
            if nargin < 2
                X = Distribution.DefaultPdfX;
                obj.isdefault = true;
            end
            switch (obj.Type)
                case DType.exp
                    lambda = obj.Params(1);
%                     assert(lambda > 0);
                    assert(lambda ~=0);
%                     if nargin < 2
%                         maxX = - 1/lambda * log( obj.PdfTol / lambda );
%                         X = linspace(0, maxX, obj.PdfGridNumber);
%                     end
                    if lambda > 0
                        domain = X>=0; 
                    else
                        domain = X<=0;
                    end
                    obj.PdfX = X;
                    obj.PdfY = domain .* (abs(lambda) * exp(-lambda * X));
                case DType.gamma
                    p = obj.Params(1);
                    b = obj.Params(2);
                    assert(p>0 && b>0);
                    obj.PdfX = X;
                    if mod(p,1) ~= 0 %|| p == 1% is not an integer?
                        fak = b^p / gamma(p);
                        obj.PdfY = fak * (X>=0) .* X .^(p-1) .* exp(-b*X);
                    else % is integer:
                        % stable version:
                        ebX = exp(-b*X);
                        Y = b * ones(size(X)) .* (X>=0) .* ebX;
                        for i=1:(p-1)
                            Y = Y .* (b/i) .* X;
                        end
                        obj.PdfY = Y;
                        % monotony check:
                        dYdX = diff(Y) ./ diff(X);
                        indp = dYdX >  eps;
                        indn = dYdX < -eps;
                        mbool = false;
                        if ~all(~indp) && ~all(~indn)
                            if find(indp, 1, 'last') > find(indn, 1, 'first')
                                if Distribution.EnableDistributionWarnings
                                    warning('gamma(%g,%g) monotony check failed', p, b);
                                end
                                mbool = true;
                            end
                        end
                        
                        % mean value check:
                        vbool = false;
                        m = p/b;
                        v = p/b^2;
                        xi = find(X>=m, 1, 'first');
                        if Y(xi) < 1 / (6*sqrt(v))
                            if Distribution.EnableDistributionWarnings 
                                warning('gamma(%g,%g) mean value check failed (Y(xmean) = %g)', p, b, Y(xi));
                            end
                            vbool = true;
                        end
                        
                        if Distribution.EnableDistributionWarnings && m < X(1) || m > X(end)
                            warning('domain too small! mean = %g, var = %g', m, v); 
                        end
                        
                        %
%                         icrit = find(ebX > 0, 1, 'last');
                        if mbool || vbool %Y(icrit) > 1e-3 % any(isnan(Y)) 
                            % recursively reduce p
%                             if Distribution.EnableDistributionWarnings 
                            if Distribution.DistributionWarningLvl > 1
                                warning('Recursively determine gamma dist. (%d --> %d + %d)', p, floor(p/2), ceil(p/2));
                            end
                            Dhelp1 = Distribution(DType.gamma, floor(p/2), b);
                            if  floor(p/2) == ceil(p/2)
                                Dhelp2 = Dhelp1;
                            else
                                Dhelp2 = Distribution(DType.gamma, ceil(p/2), b);
                            end
                            obj = Dhelp1 + Dhelp2;
                        end
                    end
                case DType.beta
                    alp = obj.Params(1);
                    bet = obj.Params(2);
                    a   = obj.Params(3);
                    b   = obj.Params(4);
                    obj.PdfX = X;
                    obj.PdfY = (X>a & X < b) .* (X-a).^(alp-1) .* (b-X).^(bet-1) / (beta(alp, bet) * (b-a)^(alp+bet-1) );
                case DType.gauss
                    mu = obj.Params(1);
                    sigmasq = obj.Params(2);
%                     if nargin < 2
%                         maxX = - 1/lambda * log( obj.PdfTol / lambda );
%                         X = linspace(0, maxX, obj.PdfGridNumber);
%                     end
                    obj.PdfX = X;
                    obj.PdfY = 1 / sqrt(2 * pi * sigmasq) * exp(- (X-mu).^2 / (2*sigmasq));
                case DType.mean
                    
                case DType.deterministic
                    obj.PdfX = X;
                    obj.PdfY = zeros(size(X));
                case DType.numeric
                    if nargin > 1
                    % Interpolate pdf on the new domain:
%                         if Distribution.EnableDistributionWarnings 
                        if Distribution.DistributionWarningLvl > 1
                            warning('Interpolation...');
                        end
%                         obj.PdfY = interp1( obj.PdfX, obj.PdfY, X, 'spline', 0);
                        obj.PdfY = interp1( obj.PdfX, obj.PdfY, X, 'linear', 0); % preserves positive values
                        obj.PdfX = X;
                    end
                case DType.uniform
                    l = obj.Params(1);
                    r = obj.Params(2);                    
                    obj.PdfX = X;
                    obj.PdfY = 1/(r-l) * (X>=l & X<r) ;
                case DType.fishertippet
                    mu = obj.Params(1);
                    sigma = obj.Params(2);
                    xi = obj.Params(3);
                    obj.PdfX = X;
                    if xi == 0
                        tx = exp(-(X-mu)/sigma);
                    else
                        tx = (1+xi*(X-mu)/sigma).^(-1/xi);
                    end
                    obj.PdfY = (X > obj.SupportLeft & X < obj.SupportRight) .* 1/sigma .* tx.^(xi+1) .* exp(-tx);
                case DType.waitingDet
                    lambda = obj.Params(1);
                    servicetime = obj.Params(2);
                    c = obj.Params(3);
                    assert(lambda > 0 && servicetime > 0 && c == 1);
                    
                    [~, PDF, p] = waitingMDc(X, lambda, servicetime, c) ;
                    obj.PdfX = X;
                    obj.PdfY = PDF;
                    obj.DetermX = 0;
                    obj.DetermY = p(1);
                case DType.invalid
%                     if Distribution.EnableDistributionWarnings
                    if Distribution.DistributionWarningLvl > 1
                        warning('invalid distribution!')
                    end
            end
            % check consistency:
            obj = checkbounds(obj);
        end
        %% --------------------------------------------------------------- %
        function obj = updateCdf(obj)
            if Distribution.DebugMode
                fprintf('DEBUG: updateCdf([%s])\n', sprintf('%d ',  obj.ID));
            end
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                res = Distribution.initDist(size(obj));
                for oi=1:numel(res)
                    res(oi) = obj(oi).updateCdf();
                end
                obj = res;
                return;
            end
            % 
            
%             warning('no sophisticated integration scheme yet');
            if obj.Type == DType.invalid 
%                 obj = Distribution(); %return invalid
                return;
            end

            obj.CdfX = obj.PdfX;
            
            indL = obj.PdfX < obj.SupportLeft;
            indR = obj.PdfX > obj.SupportRight;
            ind  = ~indL & ~indR;
            
            switch (obj.Type)
                case DType.exp
                    lambda = obj.Params(1);                    
                    if lambda > 0
                        domain = obj.CdfX>=0;                         
                        obj.CdfY = domain .* (1 - exp(- lambda * obj.CdfX));
                    else
                        domain = obj.CdfX<=0;                        
                        obj.CdfY = domain .* exp(- lambda * obj.CdfX);
                        obj.CdfY(~domain) = 1;
                    end                    
%                 case DType.gauss
%                     mu = obj.Params(1);
%                     sigmasq = obj.Params(2);
% %                     if nargin < 2
% %                         maxX = - 1/lambda * log( obj.PdfTol / lambda );
% %                         X = linspace(0, maxX, obj.PdfGridNumber);
% %                     end
%                     obj.PdfY = 1 / sqrt(2 * pi * sigmasq) * exp(- (obj.CdfX-mu).^2 / (2*sigmasq));
%                 case DType.mean
                    
%                 case DType.numeric
                case DType.beta
                    alp = obj.Params(1);
                    bet = obj.Params(2);
                    a   = obj.Params(3);
                    b   = obj.Params(4);
                    X = obj.CdfX;
                    obj.CdfY(X<a) = 0;
                    obj.CdfY(X>=b)= 1;
                    ind           = (X >= a) & (X<b);
                    X             = X(ind);
                    obj.CdfY(ind) = betainc((X-a)/(b-a), alp, bet);
                case DType.gamma
                    p = obj.Params(1);
                    b = obj.Params(2);    
                    X = obj.CdfX;
                    obj.CdfY(X<0) = 0;
                    obj.CdfY(X>=0)= gammainc(b*X(X>=0), p);
                case DType.uniform
                    l = obj.Params(1);
                    r = obj.Params(2);                    
                    X = obj.CdfX;
                    obj.CdfY = 1/(r-l) * (X>=l & X<r) .* (X-l) + (X>=r);
                case DType.fishertippet
                    mu = obj.Params(1);
                    sigma = obj.Params(2);
                    xi = obj.Params(3);
                    X = obj.CdfX;
                    if xi == 0
                        tx = exp(-(X-mu)/sigma);
                    else
                        tx = (1+xi*(X-mu)/sigma).^(-1/xi);
                    end
                    obj.CdfY = (X > obj.SupportLeft & X < obj.SupportRight) .* exp(-tx);
                    obj.CdfY(X>=obj.SupportRight) = 1;
                    
                case DType.waitingDet
                    lambda = obj.Params(1);
                    servicetime = obj.Params(2);
                    c = obj.Params(3);
                    assert(lambda > 0 && servicetime > 0 && c == 1);
                    
                    X = obj.CdfX;
                    [CDF, ~, ~] = waitingMDc(X, lambda, servicetime, c) ;
                    obj.CdfY = CDF;
                    
                case DType.invalid
%                     if Distribution.EnableDistributionWarnings 
                    if Distribution.DistributionWarningLvl > 1
                        warning('invalid distribution!')
                    end
                    
                otherwise % numeric:
%                     dx = diff(obj.CdfX);
                    % Mittelwertregel:
%                     yy = 0.5 * (obj.PdfY(1:end-1) + obj.PdfY(2:end));
%                     obj.CdfY = [0 cumsum(dx .* yy)];
                    
                    
                    obj.CdfY = zeros(size(obj.PdfX));
                    if sum(ind) > 1
                        obj.CdfY(ind)  = numericCDF(obj.PdfX(ind), obj.PdfY(ind), 'NC', Distribution.IntegrationPoints);
                    end
            end
            
            % add deterministic parts to CDF:
            if obj.Type ~= DType.waitingDet
                % (deterministic parts are already included in waitingDet)
                dx = diff(obj.CdfX(1:2));
                indadd = obj.DetermX > obj.CdfX(end);
                DX = obj.DetermX(indadd);
                if ~isempty(DX)
                    obj.CdfX = [obj.CdfX sort([DX-dx, DX])];
                    obj.CdfY = [obj.CdfY obj.CdfY(end) * ones(1, 2*sum(indadd))];
                end
                for i=1:length(obj.DetermX)
                    obj.CdfY( obj.CdfX >= obj.DetermX(i) ) = obj.CdfY( obj.CdfX >= obj.DetermX(i) ) + obj.DetermY(i);
                end
            end
            
            % normalize, if necessary:
%             normfak     = obj.CdfY(find(ind, 1, 'last'));
%             if normfak ~= 0
%                 obj.PdfY    = obj.PdfY / normfak;
%                 obj.DetermY = obj.DetermY / normfak;
%                 obj.CdfY    = obj.CdfY / normfak;
%             end
            % ---- %
            
            obj = obj.normalize();
            
            if obj.NonNormalized
                obj.CdfY(indL) = min(obj.CdfY);
                obj.CdfY(indR) = max(obj.CdfY);
            else
                obj.CdfY(indL) = 0;
                obj.CdfY(indR) = 1;
            end
            obj = checkbounds(obj);
        end
        % --------------------------------------------------------------- %
        function obj = changeXDomain(obj, newxdomain)
            if Distribution.DebugMode
                fprintf('DEBUG: changeXDomain([%s])\n', sprintf('%d ', obj.ID));
            end
            %already handled in updatePDF
            obj = updatePdf(obj, newxdomain);
            obj = updateCdf(obj);
        end
        %% --------------------------------------------------------------- %
        function obj = normalize(obj)
            if Distribution.DebugMode
                fprintf('DEBUG: normalize([%s])\n', sprintf('%d ', obj.ID));
            end
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                res = Distribution.initDist(size(obj));
                for oi=1:numel(res)
                    res(oi) = obj(oi).normalize();
                end
                obj = res;
                return;
            end
            %
            
            if obj.Type == DType.invalid
                obj = Distribution(); %return invalid
                return;
            end
            
            if obj.NonNormalized
                if Distribution.DistributionWarningLvl > 1
                    warning('NonNormalized Distribution');
                end
                return;
            end
            indL = obj.PdfX <= obj.SupportLeft;
            indR = obj.PdfX >= obj.SupportRight;
            ind  = ~indL & ~indR;
            normfak     = obj.CdfY(find(ind, 1, 'last'));
            if normfak ~= 0
                obj.PdfY    = obj.PdfY / normfak;
                obj.DetermY = obj.DetermY / normfak;
                obj.CdfY    = obj.CdfY / normfak;
            end
            obj = checkbounds(obj);
        end
        %% --------------------------------------------------------------- %
        function obj = checkbounds(obj)
            if Distribution.DebugMode
                fprintf('DEBUG: checkbounds([%s])\n', sprintf('%d ', obj.ID));
            end
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                res = Distribution.initDist(size(obj));
                for oi=1:numel(res)
                    res(oi) = obj(oi).checkbounds();
                end
                obj = res;
                return;
            end
            %
            
            % invalid:
            if obj.Type == DType.invalid 
                obj = Distribution(); %return invalid
                return
            end
            
            if obj.NonNormalized
                if Distribution.DistributionWarningLvl > 1
                    warning('NonNormalized Distribution');
                end
                return;
            end
            
            if any(obj.PdfY < 0)
                err = min(obj.PdfY);
                if Distribution.EnableDistributionWarnings  
                    if Distribution.DistributionWarningLvl == 2 || abs(err) > Distribution.DistributionWarningTol
                        warning('pdf < 0 detected. Correct error of %g', err);
                    end
                end
                obj.PdfY = max(obj.PdfY, 0);
            end
            % CDF bounds:            
            if any(obj.CdfY < 0)
                err = min(obj.CdfY);
                if Distribution.EnableDistributionWarnings 
                    if Distribution.DistributionWarningLvl == 2 || abs(err) > Distribution.DistributionWarningTol
                        warning('cdf < 0 detected. Correct error of %g', err);
                    end
                end
                obj.CdfY = max(obj.CdfY, 0);
            end
            if any(obj.CdfY > 1)
                err = max(obj.CdfY-1);
                if Distribution.EnableDistributionWarnings 
                    if Distribution.DistributionWarningLvl == 2 || abs(err) > Distribution.DistributionWarningTol
                        warning('cdf > 1 detected. Correct error of %g', err);
                    end
                end
                obj.CdfY = min(obj.CdfY, 1);
            end
        end
        %% --------------------------------------------------------------- %
        function obj = cleanDist(obj)
            if Distribution.DebugMode
                fprintf('DEBUG: cleanDist([%s])\n', sprintf('%d ', obj.ID));
            end
            % Cleans the distribution. This includes merging close values
            % of the discrete part.  
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                res = Distribution.initDist(size(obj));
                for oi=1:numel(res)
                    res(oi) = obj(oi).cleanDist();
                end
                obj =res;
                return;
            end
            %
            
            % deterministic part:
            if ~isempty(obj.DetermX) 
                % sort, if necessary:
                [newdetx, di] = sort(obj.DetermX);
                newdety = obj.DetermY(di);

                % find "clusters":
                diffs = [diff(newdetx) inf];
                % absolute:
    %             ind   = abs(diffs) < Distribution.DiscTol1;
                % relative:
                ind   = false(size(newdetx));
                irel  = abs(newdetx) > 1;
                iabs  = ~irel;
                ind(irel) = abs(diffs(irel)) ./ newdetx(irel) < Distribution.DiscTol1;
                ind(iabs) = abs(diffs(iabs)) < Distribution.DiscTol1;


                clusterlength = zeros(size(ind));
                clusterlength(1) = ind(1);
                for i=2:length(ind)
                    if ind(i)
                        clusterlength(i) = clusterlength(i-1) + 1;
                    else
                        clusterlength(i) = 0;
                    end
                end

                % merge very close discrete values:
                inddelete   = false(1,length(newdetx));
                clusterends = find( diff(clusterlength) < 0 ); % -1 is already a cluster of two values
                insdetx = zeros(size(clusterends));
                insdety = insdetx;
                for i=1:length(clusterends)
                    pos = clusterends(i);
                    iii = 1+pos-clusterlength(pos):pos+1; % we have to add 1 to the end
                    insdetx(i) = mean( newdetx( iii ) );
                    insdety(i) = sum ( newdety( iii ) );
                    inddelete(iii) = true;
                end
                newdetx(inddelete) = [];
                newdety(inddelete) = [];

                newdetx = [newdetx insdetx];
                newdety = [newdety insdety];

                [newdetx, di] = sort(newdetx);
                newdety = newdety(di);

                % remove values with extremely low probability:
                ind = newdety < Distribution.DiscTol2;

                %add them to the pdf:
                dx = diff(obj.PdfX(1:2));
                nind = sum(ind);
                nx   = length(obj.PdfX);
                if nind > 0
                    if nind * nx * 8 < 500 * (1024)^2 % size is smaller than 500MB
                        [~, mi] = min ( abs(ones(nx, 1) * newdetx(ind) - obj.PdfX.' * ones(1, nind)) );
                    else
                        % manually due to space limitations:
                        mi = nan(1, nind);
                        for i=1:nind
                            [~, mi(i)] = min ( abs(newdetx(i) - obj.PdfX) );
                        end
                    end                
                    obj.PdfY(mi) = obj.PdfY(mi) + newdety(ind)  / dx;
                end
                %            
                newdetx(ind) = [];
                newdety(ind) = [];
                obj.DetermX = newdetx;
                obj.DetermY = newdety;      
            end
            
            obj = obj.updateCdf();
            obj = obj.normalize();
        end
        %% --------------------------------------------------------------- %
        function E = expectation(obj, f)     
            if Distribution.DebugMode
                fprintf('DEBUG: expectation([%s])\n', sprintf('%d ', obj.ID));
            end
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                E = nan(size(obj));
                for oi=1:numel(obj)
                    E(oi) = obj(oi).expectation();
                end
                return;
            end
            %
            
            % invalid:
            if obj.Type == DType.invalid 
                E = NaN; %return invalid
                return;
            end
            
            if nargin < 2
                f = @(x) x; 
            end
            indL = obj.PdfX < obj.SupportLeft;
            indR = obj.PdfX > obj.SupportRight;
            ind  = ~indL & ~indR;
            
            
            Econt = 0;
            Edisc = 0;
            if obj.Type ~= DType.deterministic
                domaincheck = numericCDF(obj.PdfX(ind), obj.PdfY(ind), 'NC', Distribution.IntegrationPoints);
                checksum = domaincheck(end) + sum(obj.DetermY);
%                 if  Distribution.EnableDistributionWarnings && checksum < 1-eps(1)
                if  Distribution.EnableDistributionWarnings && checksum < 1-Distribution.DistributionWarningTol
                    warning('Numeric Domain to small for Expectation calculation! (Int = %g)', checksum);
                end
                Econt = numericCDF(obj.PdfX(ind), f(obj.PdfX(ind)) .* obj.PdfY(ind), 'NC', Distribution.IntegrationPoints);
                Econt = Econt(end);
            end
            if ~isempty(obj.DetermX)
                Edisc = f(obj.DetermX) * obj.DetermY';
            end
            
            E = Econt + Edisc;            
        end
        %% --------------------------------------------------------------- %
        function [obj, newcdfx, newcdfy] = functionofX(obj, g, invg, dinvg)
            if Distribution.DebugMode
                fprintf('DEBUG: functionofX([%s], ...)\n', sprintf('%d ', obj.ID));
            end
            % gets the distribution of g(X)
            % if the inverse invg and derivative dinvg of the inverse is also given, an analytic solution is available
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                res = Distribution.initDist(size(obj));
                newcdfx = cell(size(obj));
                newcdfy = cell(size(obj));
                for oi=1:numel(res)
                    if nargin < 3
                        [res(oi), newcdfx{oi}, newcdfy{oi}] = obj(oi).functionofX(g);
                    else
                        [res(oi), newcdfx{oi}, newcdfy{oi}] = obj(oi).functionofX(g, invg, dinvg);
                    end
                end
                obj = res;
                return;
            end
            %
            
            if obj.Type == DType.invalid 
                obj = Distribution(); %return invalid
                newcdfx = Distribution.DefaultPdfX;
                newcdfy = newcdfx*0;
                return;
            end
            
            if nargin < 3
                analytic = false;
            else
                analytic = true;
            end
            
            if analytic
                % analytic solution:
                newpdfx = obj.PdfX;
                newpdfy = obj.interpPDF( invg(newpdfx) ) .* dinvg(newpdfx);
                newdiscrx = g(obj.DetermX);
                newdiscry = obj.DetermY;
                obj = Distribution(DType.numeric, newpdfx, newpdfy, newdiscrx, newdiscry);
                %
                newcdfx = newpdfx;
                newcdfy = obj.interpCDF( invg(newpdfx) );
            else
                % numeric approximation:
                [disc, cont, probs] = getDiscreteAndContinuousPart(obj);
                % continous part:
                newcdfx = cont.CdfX;
                dx = newcdfx(2) - newcdfx(1);
                % enlarge the domain to make sure that we have enough points:
                domainlength = newcdfx(end) - newcdfx(1);
                enlargedxdomain = linspace(newcdfx(1)-domainlength, newcdfx(end)+domainlength, 3*length(newcdfx));
                gvalues = g(enlargedxdomain);
%                 [u, ia, ic] = unique(gvalues);
%                 ginvvalues = interp1(gvalues, newcdfx, newcdfx);
                [ginvvalues, jx] = uniqueinterp(gvalues, enlargedxdomain, newcdfx, true, -inf, +inf);
                newcdfy = cont.interpCDF( ginvvalues );
                jy = 0 * jx;
                for i=1:length(jx)
                    % treat the jumps in g^-1:
                    [~, mi] = min(abs(newcdfx - jx(i)));
                    if mi>1
                        jump = newcdfy(mi) - newcdfy(mi-1);
                        newcdfy(mi:end) = newcdfy(mi:end) - jump;
                        jy(i) = jump;
                    end
                end
                newpdfx = newcdfx;
%                 newpdfy = [0, newcdfy(3:end) - newcdfy(1:end-2), 0] / (2*dx); % symmetric difference quotient
                newpdfy = [newcdfy(2:end) - newcdfy(1:end-1), 0] / (dx); % symmetric difference quotient from right direction
                newcont = Distribution(DType.numeric, newpdfx, newpdfy, jx, jy);
                % discrete part:
                newdiscrx = g(disc.DetermX);
                newdiscry = disc.DetermY; 
                if isempty(disc.DetermX)
                    newdisc = Distribution();
                else               
                    newdisc = Distribution(DType.deterministic, newdiscrx, newdiscry);
                end
                % put it together:
                obj = probs * [newdisc, newcont];
                
                % update CDF parts (for single output):
                for i=1:length(newdiscrx)
                    newcdfy( newcdfx >= newdiscrx(i) ) = newcdfy( newcdfx >= newdiscrx(i) ) + newdiscry(i);
                end
            end
        end
        %% --------------------------------------------------------------- %
        function obj = functionofXandY(X,Y,g)
            % gets the distribution of g(X,Y) for independent RV X and Y.
            % g should be implemented in a way that it can operate on elementwise on matrices
            if Distribution.DebugMode
                fprintf('DEBUG: functionofXandY([%s], [%s], g)\n', sprintf('%d ', X.ID), sprintf('%d ', Y.ID));
            end
            
            % no array functionality so far
            
            if X.Type == DType.invalid || Y.Type == DType.invalid
                obj = Distribution(); %return invalid
                return;
            end
            
            [Xdisc, Xcont, Xprobs] = getDiscreteAndContinuousPart(X);
            [Ydisc, Ycont, Yprobs] = getDiscreteAndContinuousPart(Y);
            
            W = [ Xprobs(1)*Yprobs(1), Xprobs(1)*Yprobs(2), Xprobs(2)*Yprobs(1), Xprobs(2)*Yprobs(2)];            
            nx = length(Xdisc.DetermX);
            ny = length(Ydisc.DetermY);
            % discrete Part:           
            if W(1) > 0
                XX = ones(ny, 1) * Xdisc.DetermX;
                YY = (Ydisc.DetermX).' * ones(1,nx);
                GG = g(XX,YY);
                PXX = ones(ny, 1) * Xdisc.DetermY;
                PYY = (Ydisc.DetermY).' * ones(1,nx);
                PGG= PXX.*PYY;
                newdiscx = GG(:); % double entries will be removed by constructor
                newdiscy = PGG(:); 
                Ddd = Distribution(DType.deterministic, newdiscx, newdiscy);
            else
                Ddd = Distribution();
            end
            % hybrid Parts:
            if W(2) > 0
                w = Xdisc.DetermY;
                Dx(nx) = Distribution();
                for i=1:nx
                    g2 = @(y) g(Xdisc.DetermX(i), y);
                    Dx(i) = functionofX(Ycont, g2);
                end
                Ddc = w * Dx;
            else
                Ddc = Distribution();
            end
            if W(3) > 0
                w = Ydisc.DetermY;
                Dy(ny) = Distribution();
                for i=1:ny
                    g2 = @(x) g(x, Ydisc.DetermX(i));
                    Dy(i) = functionofX(Xcont, g2);
                end
                Dcd = w * Dy;
            else
                Dcd = Distribution();                
            end
            % continuous Parts:
            if W(4) > 0
                newcdfx = Xcont.CdfX;
                dx = newcdfx(2) - newcdfx(1);
                newcdfy = 0 * newcdfx;
                ind = find( Xcont.PdfX > 1e-9);
%                 for i=ind
                for j=1:length(ind)
                    if mod(j, 100)==0
                        fprintf( '%d / %d\n', j, length(ind));
                    end
                    i = ind(j);
                    g2 = @(y) g(Xcont.PdfX(i), y); % x fix
                    [~, FYx, FYy] = Ycont.functionofX(g2);
                    newcdfy = newcdfy + Xcont.PdfY(i) * FYy;
                end                 
                newcdfy = newcdfy * dx;
                newpdfx = newcdfx;
                newpdfy = [0, newcdfy(3:end) - newcdfy(1:end-2), 0] / (2*dx); % symmetric difference quotient
%                 newpdfy = [newcdfy(2:end) - newcdfy(1:end-1), 0] / (dx); % symmetric difference quotient from right direction
                Dcc = Distribution(DType.numeric, newpdfx, newpdfy);
            else
                Dcc = Distribution();
            end
            obj = W * [Ddd, Ddc, Dcd, Dcc];
        end
        %% --------------------------------------------------------------- %
        function obj = minus(D1, D2)
            % Distribution of the difference of both (independent) random variables
            % given by their distributions D1 and D2.
            % 
            % just reuse + and handle special cases.
            if Distribution.DebugMode
                fprintf('DEBUG: minus([%s], [%s])\n', sprintf('%d', D1.ID), sprintf('%d ', D2.ID));
            end
            
            % Handling of arrays of Distributions:
            if numel(D1) > 1
                obj = Distribution.initDist(size(D1));
                if numel(D2) == numel(D1)
                    selector = @(i) i;
                elseif numel(D2) == 1
                    selector = @(i) 1;
                else
                    error('In D1-D2 number of elements must be the same or scalar!');
                end
                %    
                for oi=1:numel(obj)                    
                    obj(oi) = D1(oi)-D2(selector(oi));
                end
                return;
            elseif numel(D2) > 1 % only the second one is an array:
                obj = Distribution.initDist(size(D2));
                for oi=1:numel(obj)                    
                    obj(oi) = D1(1)-D2(oi);
                end
                return;
            end
            %
            
            % subtract a deterministic (scalar) value:
            if isa(D2, 'double')
                obj = D1 + (-D2);
                return
            end
            
            if D1.Type == DType.invalid || D2.Type == DType.invalid
                obj = Distribution(); %return invalid
                return;
            end
            
            % standard case:
            obj = D1 + D2.multiply(-1);            
        end
        %% --------------------------------------------------------------- %
        function obj = uminus(D1)
            % Distribution of a negated random variable
            % given by its distribution D1.
            % 
            % just reuse the multiply function.
            if Distribution.DebugMode
                fprintf('DEBUG: uminus([%s])\n', sprintf('%d ', D1.ID));
            end
            
            obj = D1.multiply(-1);
        end
        %% --------------------------------------------------------------- %
        function obj = plus(D1, D2)
            % Distribution of the sum of both (independent) random variables
            % given by their distributions D1 and D2.
            % = convolution of the pdf's
%             warning('plus not yet implemented');
%             obj = Distribution(DType.invalid);
            if Distribution.DebugMode
                fprintf('DEBUG: plus([%s], [%s])\n', sprintf('%d', D1.ID), sprintf('%d ', D2.ID));
            end                       
            

            % add a deterministic (scalar) value:
            if isa(D1, 'double')
                obj = D2.shift(D1);
                return
            end
            if isa(D2, 'double')
                obj = D1.shift(D2);
                return
            end
            
            % Handling of arrays of Distributions:
            if numel(D1) > 1
                obj = Distribution.initDist(size(D1));
                if numel(D2) == numel(D1)
                    selector = @(i) i;
                elseif numel(D2) == 1
                    selector = @(i) 1;
                else
                    error('In D1+D2 number of elements must be the same or scalar!');
                end
                %    
                for oi=1:numel(obj)                    
                    obj(oi) = D1(oi)+D2(selector(oi));
                end
                return;
            elseif numel(D2) > 1 % only the second one is an array:
                obj = Distribution.initDist(size(D2));
                for oi=1:numel(obj)                    
                    obj(oi) = D1(1)+D2(oi);
                end
                return;
            end
            %
            
            % invalid:
            if D1.Type == DType.invalid || D2.Type == DType.invalid
                obj = Distribution(); %return invalid
                return
            end
            
            %
            resultNonNormalized = D1.NonNormalized || D2.NonNormalized;            
            
            % deterministic value:
            if D1.Type == DType.deterministic
                if length(D1.DetermX) == 1
                    obj = D2.shift(D1.DetermX);
                    return
                end
            end
            if D2.Type == DType.deterministic
                if length(D2.DetermX) == 1
                    obj = D1.shift(D2.DetermX);
                    return
                end
            end
            
            % in case nothing else is possible:
%             warning('ensure that domains are the same and equally distributed');
            
            done = false;
%             if D1.Type == DType.gauss && D2.Type == DType.gauss
%                 mu = D1.Params(1) + D2.Params(1);
%                 sigmasq = D1.Params(2) + D2.Params(2);
%                 obj = Distribution(DType.gauss, mu, sigmasq);
%                 done = true;
%             end

            sl = D1.SupportLeft  + D2.SupportLeft;
            sr = D1.SupportRight + D2.SupportRight;
            
            % two deterministic:
            if (D1.Type == DType.deterministic) && (D2.Type == DType.deterministic)
                % values:
                determX = D1.DetermX' * ones(size(D2.DetermX)) + ones(size(D1.DetermX))' * D2.DetermX;
                determX = determX(:)';
                %probs:
                determY = D1.DetermY' * D2.DetermY; % all combinations
                determY = determY(:)';
                obj = Distribution(DType.deterministic, determX, determY, resultNonNormalized);
%                 done = true;
                return
            end
            
            % two gaussians:
            if (D1.Type == DType.gauss) && (D2.Type == DType.gauss)
                mu = D1.Params(1)+D2.Params(1);
                sigma2 = D1.Params(2) + D2.Params(2);
                obj = Distribution(DType.gauss, mu, sigma2);
%                 done = true;
                return
            end
            
            % two exponentials / gammas
            if (D1.Type == DType.exp || D1.Type == DType.gamma) && (D2.Type == DType.exp || D2.Type == DType.gamma)
                if (D1.Type == DType.exp)                    
                    p1 = 1;
                    b1 = D1.Params(1);
                else
                    p1 = D1.Params(1);
                    b1 = D1.Params(2);
                end
                if (D2.Type == DType.exp)                    
                    p2 = 1;
                    b2 = D2.Params(1);
                else
                    p2 = D2.Params(1);
                    b2 = D2.Params(2);
                end
                
                if b1 == b2 && b1 > 0 % only then simplification is possible
                    obj = Distribution(DType.gamma, p1+p2, b1);
    %                 done = true;
                    return
                end
            end
            
            % two uniform
            if (D1.Type == DType.uniform) && (D2.Type == DType.uniform)
                if Distribution.EnableDistributionWarnings && Distribution.DistributionWarningLvl > 1
                    warning('t.b.i.: Sum of two uniform distributions');
                end
            end
                        
            if ~done
                detn1 = length(D1.DetermX);
                detn2 = length(D2.DetermX);
%                 p1 = sum(D1.DetermY);
%                 p2 = sum(D2.DetermY);
                W  = [];
                D  = [];
                % deterministicpart is included in the hybrid part:
%                 if  ~isempty(D1.DetermX) && ~isempty(D2.DetermX)
%                     % values:
%                     determX = D1.DetermX' * ones(size(D2.DetermX)) + ones(size(D1.DetermX))' * D2.DetermX;
%                     determX = determX(:);
%                     %probs:
%                     determY = D1.DetermY' * D2.DetermY; % all combinations
%                     determY = determY(:);
%                     W = [W 1];
%                     D = [D Distribution(DType.deterministic, determX, determY)];
%                 end
                % hybrid part:
                if ~isempty(D1.DetermX)
                    % old single version:
                    Dtmp(detn1) = Distribution(DType.invalid);
                    for i=1:detn1
                        Dtmp(i) = D2.shift(D1.DetermX(i)); %shift deterministic and numeric part
                    end
%                     W = [W 1];
%                     D = [D (D1.DetermY * Dtmp)];
                    W = [W D1.DetermY];
                    D = [D Dtmp];
                    % new batch version:
%                     Dtmp = D2.shift(D1.DetermX, D1.DetermY);
%                     W = [W sum(D1.DetermY)];
%                     D = [D Dtmp];
                    clear Dtmp;
                end
                if ~isempty(D2.DetermX) && sum(D1.DetermY) < 1 % && ~(D1.Type == DType.deterministic)
                    % take only the continuous part. (discrete part already
                    % handled)
                    D1cont = D1;                    
                    D1cont.DetermX = []; 
                    D1cont.DetermY = [];
                    D1cont.PdfY = 1 / (1-sum(D1.DetermY)) * D1cont.PdfY;
                    % old single version:
                    Dtmp(detn2) = Distribution(DType.invalid);
                    for i=1:detn2
                        Dtmp(i) = D1cont.shift(D2.DetermX(i)); %shift deterministic and numeric part
                    end
%                     W = [W 1];
%                     D = [D (D2.DetermY * Dtmp)];
                    W = [W (1-sum(D1.DetermY)) * D2.DetermY];
                    D = [D Dtmp];
                    % new batch version:
%                     Dtmp = D1cont.shift(D2.DetermX, D2.DetermY);
%                     W = [W (1-sum(D1.DetermY)) * sum(D2.DetermY)];
%                     D = [D Dtmp];
                    clear Dtmp;
                end
                % continous part: (both need a continuous part to let the
                %                 following have a sense.)
                if ~(D1.Type == DType.deterministic) && ~(D2.Type == DType.deterministic)
                    % option 1:
                    w1 = 1-sum(D1.DetermY);
                    w2 = 1-sum(D2.DetermY);
                    if w1 == 0 || w2 == 0
                        X = D1.PdfX;
                        Y = 0 * X;
                    else
                        [X, Y] = convolve(D1.PdfX, D1.PdfY / w1, D2.PdfX, D2.PdfY / w2);
                    end
%                 obj = Distribution(DType.numeric, X, Y, determX, determY);
                    W = [W 1-sum(W)];
                    % option 2:
%                     [X, Y] = convolve(D1.PdfX, D1.PdfY, D2.PdfX, D2.PdfY);
% %                 obj = Distribution(DType.numeric, X, Y, determX, determY);
%                     W = [W 1];
                    D = [D Distribution(DType.numeric, X, Y, [], [], resultNonNormalized)];
                end
                obj = W * D;
                
                % it is very likely that the discrete part will be crowded
                % afterwards
                obj = obj.cleanDist();
                
                % update support:
                obj.SupportLeft = sl;
                obj.SupportRight = sr;
                
                % part of cleanDist():
%                 obj = obj.updateCdf;                
%                 obj = obj.normalize();
            end
        end
        %% --------------------------------------------------------------- %
        function obj = or(D1, D2)
            % Distribution of the alternative of both Distributions
            % no weights specified --> 0.5 for each
%             warning('or not yet tested');
%             obj = Distribution(DType.invalid);
            if Distribution.DebugMode
                fprintf('DEBUG: or(%d, %d)\n', D1.ID, D2.ID);
            end
            
            if numel(D1) > 1 || numel(D2) > 1
                warning('mtimes: No Array Handling implemented yet');
            end

%             support = [min(D1.Support(1), D2.Support(1)) max(D1.Support(2), D2.Support(2))];
            sl = min([D1.SupportLeft D2.SupportLeft]);
            sr = max([D1.SupportRight D2.SupportRight]);

            % in case nothing else is possible:
%             if D1.isdefault && D2.isdefault
%                 X = Distribution.DefaultPdfX;
%                 Y1= D1.PdfY;
%                 Y2= D2.PdfY;                
%             elseif isequal(D1.PdfX, D2.PdfX)
%                 X = D1.PdfX;
%                 Y1= D1.PdfY;
%                 Y2= D2.PdfY;
%             else
%                 warning('Pdf domains differ. interpolation necessary');
%                 X = sort(unique([D1.PdfX D2.PdfX]));
%                 Y1 = interp1(D1.PdfX, D1.PdfY, X, 'spline', 0);
%                 Y2 = interp1(D2.PdfX, D2.PdfY, X, 'spline', 0);
%             end
%             Y = 0.5 * (Y1 + Y2);
%             obj = Distribution(DType.numeric, X, Y);            
            %replace by mtimes:
            obj = [0.5 0.5] * [D1 D2];
            % update support:
            obj.SupportLeft = sl;
            obj.SupportRight = sr;
            obj = obj.updateCdf;
        end
        %% --------------------------------------------------------------- %
        function obj = mtimes(weights, D)
            % vector of weights times vector of distributions
%             obj = Distribution(DType.invalid);
            if Distribution.DebugMode
                fprintf('DEBUG: mtimes(%s)\n', fprintf('%d ', [D.ID]) );
            end
            
            if isa(D, 'double') 
                % try swapping:
                tmp = D;
                D = weights;
                weights = tmp;
            end
            assert(isa(weights, 'double') && isa(D, 'Distribution'));
            
            if size(D,1) > 1 && size(D,2) > 1
                warning('mtimes: No Array Handling implemented yet');
            end
                        
            resultNonNormalized = any( [ D.NonNormalized ] );
            
            weightserr = abs(sum(weights)-1);
            if ~resultNonNormalized && Distribution.EnableDistributionWarnings && weightserr > Distribution.DistributionWarningTol
                warning('weights do not sum up to 1. Error: %g\n', weightserr); 
            end
            
            
            if numel(weights) == 1
                obj = D;
                return
            end
            
            % remove Distributions with zero probability
            iszero = weights == 0;
            weights(iszero)=[];
            D(iszero) = [];
            
%             support = [min(D1.Support(1), D2.Support(1)) max(D1.Support(2), D2.Support(2))];
            sl = min([D.SupportLeft]);
            sr = max([D.SupportRight]);
            
            alldet = all([D.Type] == DType.deterministic);
            
            % if any of the Distributions with nonzero probability is
            % invalid the result will be invalid
            if any( [D.Type] == DType.invalid )
                obj = Distribution(); %return invalid
                return;
            end

            % in case nothing else is possible:
            
            %deterministic part:
            determX = [];
            determY = [];
            for i = 1:length(D)
                dx = D(i).DetermX;
                dy = D(i).DetermY;
                determX = [determX dx];
                determY = [determY weights(i)*dy];
            end
            
            if alldet
                obj = Distribution(DType.deterministic, determX, determY, resultNonNormalized);
            else
                % continuous part
                if all([D.isdefault])
                    X = Distribution.DefaultPdfX;
                    Ys = reshape([D.PdfY], [], length(D))';
                else
                    %check for equal domains
                    eqx = true;
                    for i=1:length(D)-1
                        if ~isequal(D(i).PdfX, D(i+1).PdfX)
                            eqx = false;
                            break
                        end
                    end

                    if eqx
                        X = D(1).PdfX;
                        Ys = reshape([D.PdfY], [], length(D))';
                    else
                        if Distribution.EnableDistributionWarnings && Distribution.DistributionWarningLvl > 1
                            warning('Pdf domains differ. interpolation necessary');
                        end
                        X = sort(unique([D.PdfX]));
                        Ys = zeros(length(D),length(X));
                        for i = 1:length(D)
                            Ys(i,:) = interp1(D(i).PdfX, D(i).PdfY, X, 'spline', 0)';
                        end                    
                    end
                end
                Y = weights * Ys;
                % put all together:
                if isempty(determX) 
                    obj = Distribution(DType.numeric, X, Y, [], [], resultNonNormalized);                
                else
                    obj = Distribution(DType.numeric, X, Y, determX, determY, resultNonNormalized);                
                end
            end
            % update support:
            obj.SupportLeft = sl;
            obj.SupportRight = sr;
            obj = obj.updateCdf;
        end
        %% --------------------------------------------------------------- %
        function prob = lt(D1, D2)
            % Returns the probability of one RV with distribution D1 is
            % less than an independent RV with distribution D2
            % X1 < X2
            D = D1 - D2;            
            prob = D.interpCDF(0);            
            % remove the probability of the discrete part at 0:
            prob = prob - sum( D.DetermY(D.DetermX == 0) );            
        end
        %% --------------------------------------------------------------- %
        function prob = le(D1, D2)
            % Returns the probability of one RV with distribution D1 is
            % less or equal than an independent RV with distribution D2
            % X1 <= X2
            D = D1 - D2;            
            prob = D.interpCDF(0);            
        end
        %% --------------------------------------------------------------- %
        function prob = gt(D1, D2)
            % Returns the probability of one RV with distribution D1 is
            % greater than an independent RV with distribution D2
            % X1 > X2 <=> X2 < X1
            prob = D2 < D1;
        end
        %% --------------------------------------------------------------- %
        function prob = ge(D1, D2)
            % Returns the probability of one RV with distribution D1 is
            % greater or equal than an independent RV with distribution D2
            % X1 >= X2 <=> X2 <= X1          
            prob = D2 <= D1;
        end
        %% --------------------------------------------------------------- %
        function obj = shift(obj, dx, w)
            % add a scalar to a random variable
            % (shifts the pdf/cdf to the right (dx>0) or to the left (dx < 0)
            if Distribution.DebugMode
                fprintf('DEBUG: shift([%s], [%s],...)\n', sprintf('%d ', obj.ID), sprintf('%g ', dx));
            end
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                res = Distribution.initDist(size(obj));
                for oi=1:numel(obj)
                    if nargin < 3
                        res(oi) = obj(oi).shift(dx);
                    else
                        res(oi) = obj(oi).shift(dx,w);
                    end
                end
                obj = res;
                return;
            end
            %
            
            if obj.Type == DType.invalid 
                obj = Distribution(); %return invalid
                return;
            end
            
            if nargin < 3 || isempty(w)
                w = ones(size(dx)) / numel(dx);
            end
            
            % Variante 1:
%             obj.PdfX = obj.PdfX + dx;
%             obj.isdefault = false;
            
            ldx = length(dx);
            ldisc = length(obj.DetermX);
            lpdfx = length(obj.PdfX);
            
            if ldx == 1
                % scalar shift
                
                % special cases:
                if dx == 0
                    %nothing to do
                    return 
                end

                % update support:
                obj.SupportLeft  = obj.SupportLeft  + dx;
                obj.SupportRight = obj.SupportRight + dx;

                if (obj.Type == DType.uniform)
                    obj.Params = obj.Params + dx;
                    obj = obj.updatePdf();
                    obj = obj.updateCdf();
                    return
                end

                if (obj.Type == DType.gauss)
                    % shift mean 
                    obj.Params(1) = obj.Params(1) + dx;
                    obj = obj.updatePdf();
                    obj = obj.updateCdf();
                    return
                end

                % standard case:

                % update deterministic parts:
                obj.DetermX = obj.DetermX + dx;

                % Variante 2: (keep original domain)
                if ~(obj.Type == DType.deterministic)
                    % special cases will loose their special properties:
                    obj.Type = DType.numeric;
                    obj.PdfY = interp1(obj.PdfX, obj.PdfY, obj.PdfX - dx, 'linear', 0);
                end

                % Update CDF:
                obj = obj.updateCdf();
            else
                % batch shift: (hopefully faster for the plus operation)
                
                if size(dx,1) > size(dx,2)
                    % row vector:
                    dx = dx .';
                end
                if size(w,1) > size(w,2)
                    % row vector:
                    w = w .';
                end

                % update support:
                obj.SupportLeft  = obj.SupportLeft  + min(dx);
                obj.SupportRight = obj.SupportRight + max(dx);

                % standard case:

                % update deterministic parts:
%                 obj.DetermX = obj.DetermX + dx;
                if ldisc > 0
                    obj.DetermX = ones(ldx, 1) * obj.DetermX  + dx.' * ones(1, ldisc) ;
                    obj.DetermY = w .' * obj.DetermY;
                end
                
                % linearize:
                obj.DetermX = obj.DetermX(:) .';
                obj.DetermY = obj.DetermY(:) .';

                % Variante 2: (keep original domain)
                if ~(obj.Type == DType.deterministic)
                    % special cases will loose their special properties:
                    obj.Type = DType.numeric;
                    if ldx * lpdfx < 65536000 % 500 MB
                        XXX = ones(ldx, 1) * obj.PdfX -  dx.' * ones(1, lpdfx); % ldx X lpdfx
                        YYY = interp1(obj.PdfX, obj.PdfY, XXX, 'linear', 0);
                        obj.PdfY = w * YYY; % ldx X 1 * ldx X lpdfx
                    else
                        if Distribution.DistributionWarningLvl > 1
                            warning('very large shift operation');
                        end
                        YYY = obj.PdfY * 0;
                        for i=1:length(dx)
                            YYY = YYY + w(i) * interp1(obj.PdfX, obj.PdfY, obj.PdfX - dx, 'linear', 0);
                        end
                        obj.PdfY = YYY;
                    end
                end

                % Update CDF:                
%                 obj = obj.updateCdf();
                obj = obj.cleanDist();
                
            end
        end
        %% --------------------------------------------------------------- %
        function obj = multiply(obj, multiplier)
            % Multiply a random variable with a scalar (positive or negative)
            % if multiple Distributions are given, then the result will be
            % either:
            %   D_i * multiplier        if numel(multiplier) == 1
            % or
            %   D_i * multiplier_i      if numel(D) == numel(multiplier)   (elementwise)
            %
            
            if Distribution.DebugMode
                fprintf('DEBUG: multiply([%s], %g)\n', sprintf('%d ', obj.ID), multiplier);
            end
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                res = Distribution.initDist(size(obj));
                for oi=1:numel(res)
                    if numel(multiplier) == 1
                        res(oi) = obj(oi).multiply(multiplier);
                    elseif numel(multiplier) == numel(res)
                        res(oi) = obj(oi).multiply(multiplier(oi));
                    else
                        error('Distribution.multiply: sizes must be consistent (either numel(obj) == numel(multiplier) or numel(multiplier) == 1)');
                    end
                end
                obj = res;
                return;
            end
            %
            
            if obj.Type == DType.invalid 
                obj = Distribution(); %return invalid
                return;
            end
            
            % special cases:
            if multiplier == 0
                %becomes deterministic --> zero for sure
                obj = Distribution(DType.deterministic, 0);
                return 
            end
            if multiplier == 1
                %nothing to do
                return 
            end
                        
            % update support:
            obj.SupportLeft  = obj.SupportLeft  * multiplier;
            obj.SupportRight = obj.SupportRight * multiplier;
            
            if multiplier < 0
                % swap support boundaries in this case
                tmp = obj.SupportLeft;
                obj.SupportLeft = obj.SupportRight;
                obj.SupportRight = tmp;
            end
            
            if (obj.Type == DType.uniform)
                obj.Params = obj.Params * multiplier;
                obj.Params = sort(obj.Params); % if multiplier < 0
                obj = obj.updatePdf();
                obj = obj.updateCdf();
                return
            end
            
            if (obj.Type == DType.gauss)
                % shift mean 
                obj.Params(1) = obj.Params(1) * multiplier;
                obj.Params(2) = obj.Params(2) * multiplier^2;
                obj = obj.updatePdf();
                obj = obj.updateCdf();
                return
            end
            
            if (obj.Type == DType.exp)
                % shift mean 
                obj.Params(1) = obj.Params(1) / multiplier;
                obj = obj.updatePdf();
                obj = obj.updateCdf();
                return
            end
            
            % standard case:
                        
            % update deterministic parts:
            obj.DetermX = obj.DetermX * multiplier;
            
            % Variante 2: (keep original domain)
            if ~(obj.Type == DType.deterministic)
                % special cases will loose their special properties:
                obj.Type = DType.numeric;
                obj.PdfY = interp1(obj.PdfX, obj.PdfY / abs(multiplier), (obj.PdfX / multiplier), 'linear', 0); % linear is more stable at discontinuities and avoids over oscillations.
                dx = obj.PdfX(2) - obj.PdfX(1);
                indicator = max(obj.PdfY([1, end])) / dx; % to see the consequences for the cdf;
                if Distribution.DistributionWarningLvl > 0 && indicator > Distribution.DistributionWarningTol 
                    warning('Multiply probably exceeds the domain (Pdf starts/ends at y=%g * dx)', indicator)
                end
            end
                                    
            % Update CDF:
            obj = obj.updateCdf();
        end
        %% --------------------------------------------------------------- %
        function obj = multiplyPdfs(D1, D2)
            if Distribution.DebugMode
                fprintf('DEBUG: multiplyPdfs([%s], [%s])\n', sprintf('%d', D1.ID), sprintf('%d ', D2.ID));
            end
            
            % Handling of arrays of Distributions:
            if numel(D1) > 1
                obj = Distribution.initDist(size(D1));
                if numel(D2) == numel(D1)
                    selector = @(i) i;
                elseif numel(D2) == 1
                    selector = @(i) 1;
                else
                    error('In D1-D2 number of elements must be the same or scalar!');
                end
                %    
                for oi=1:numel(obj)                    
                    obj(oi) = D1(oi).multiplyPdfs(D2(selector(oi)));
                end
                return;
            elseif numel(D2) > 1 % only the second one is an array:
                obj = Distribution.initDist(size(D2));
                for oi=1:numel(obj)                    
                    obj(oi) = D1(1).multiplyPdfs(D2(oi));
                end
                return;
            end
            %
            
            if D1.Type == DType.invalid || D2.Type == DType.invalid
                obj = Distribution(); %return invalid
                return;
            end
            
            [Ddisc1, Dcont1, probs1] = getDiscreteAndContinuousPart(D1);
            [Ddisc2, Dcont2, probs2] = getDiscreteAndContinuousPart(D2);
            
            % new domain:
            newDiscrX = unique([Ddisc1.DetermX, Ddisc2.DetermX]);
            newDiscrY = newDiscrX * 0;
            newPdfX   = Dcont1.PdfX;
            
            for i=1:length(newDiscrX)
                detx = newDiscrX(i);
                i1 = find(Ddisc1.DetermX == detx, 1, 'first');
                i2 = find(Ddisc2.DetermX == detx, 1, 'first');
                % discrete parts:
                if ~isempty(i1) && ~isempty(i2)
                    % both have this entry:
                    newDiscrY(i) = newDiscrY(i) + probs1(1) * Ddisc1.DetermY(i1) *  probs2(1) * Ddisc2.DetermY(i1);
                end
                % hybrid parts:
                if ~isempty(i1)
                    newDiscrY(i) = newDiscrY(i) + probs1(1) * Ddisc1.DetermY(i1) *  probs2(2) * Dcont2.interpPDF(detx);
                end
                if ~isempty(i2)
                    newDiscrY(i) = newDiscrY(i) + probs2(1) * Ddisc2.DetermY(i2) *  probs1(2) * Dcont1.interpPDF(detx);
                end
            end                                
            % remove zeros:
            ind = newDiscrY==0;
            newDiscrX(ind) = [];
            newDiscrY(ind) = [];
            
            % continous parts:
            newPdfY = (probs1(2) * Dcont1.PdfY) .* (probs2(2) * Dcont2.PdfY);
            
            % create and return a new (non-normalized) Distribution Object
            obj = Distribution(DType.numeric, newPdfX, newPdfY, newDiscrX, newDiscrY, true);
        end
        %% --------------------------------------------------------------- %
        function obj = maxDist(D1, D2)
            % Distribution of max(D1, D2) for independent distributions D1, D2
            if nargin < 2 % D1 hopefully contains a vector in this case
                if length(D1) == 1
                    obj = D1; 
                    return
                end                
                ihalf = ceil(length(D1) / 2);
                D2 = D1(ihalf+1:end);
                D1 = D1(1:ihalf);
            end            
            
            if Distribution.DebugMode
                fprintf('DEBUG: maxDist([%s], [%s])\n', sprintf('%d ', D1.ID), sprintf('%d ', D2.ID));
            end
            
            if length(D1) > 1
                %recursively:
                D1 = maxDist(D1);
            end
            if length(D2) > 1
                %recursively:
                D2 = maxDist(D2);
            end
            
            if isa(D1, 'double')
                obj = pullBelowToX(D2, D1);
                return
            end
            if isa(D2, 'double')
                obj = pullBelowToX(D1, D2);
                return
            end
            
            if D1.Type == DType.invalid || D2.Type == DType.invalid
                obj = Distribution(); %return invalid
                return;
            end
            
            done = false;
            
            resultNonNormalized = D1.NonNormalized || D2.NonNormalized;
            
            % special cases:
            
            % general case:
            % support:
            
            sl = max(D1.SupportLeft, D2.SupportLeft);
            sr = max(D1.SupportRight, D2.SupportRight);
            
            if ~done
                [Ddisc1, Dcont1, probs1] = getDiscreteAndContinuousPart(D1);
                [Ddisc2, Dcont2, probs2] = getDiscreteAndContinuousPart(D2);            
                
                discreten1 = length(D1.DetermX);
                discreten2 = length(D2.DetermX);
                discrprob1 = probs1(1);
                discrprob2 = probs2(1);
                W = [discrprob1 * discrprob2, discrprob1 * (1-discrprob2), (1-discrprob1)*discrprob2, (1-discrprob1)*(1-discrprob2)];
                D(4) = Distribution();
                % discrete part:
                if W(1) > 0
                    newDetX = sort(unique([Ddisc1.DetermX, Ddisc2.DetermX]));
%                     newDetY = newDetX * 0;
                    DetX1    = [-Inf, Ddisc1.DetermX];
                    cumDetY1 = [0,    cumsum(Ddisc1.DetermY)];
                    DetX2    = [-Inf, Ddisc2.DetermX];
                    cumDetY2 = [0,    cumsum(Ddisc2.DetermY)];
                    cumMaxXY = newDetX * 0;
                    for i=1:length(cumMaxXY)
                        i1 = find(DetX1 <= newDetX(i), 1, 'last');
                        i2 = find(DetX2 <= newDetX(i), 1, 'last');
                        cumMaxXY(i) = cumDetY1(i1) * cumDetY2(i2);
                    end
                    newDetY = diff([0 cumMaxXY]);
                    D(1) = Distribution(DType.deterministic, newDetX, newDetY, resultNonNormalized);
                else
                    D(1) = Distribution(DType.deterministic, 0,1);
                end
                % hybrid part:
                if W(2) > 0 % discr from D1, cont from D2
                    wtmp = Ddisc1.DetermY;
                    dtmp(discreten1) = Distribution();
                    for i=1:discreten1
                        dtmp(i) = pullBelowToX(Dcont2, Ddisc1.DetermX(i));
                    end
                    D(2) = wtmp * dtmp;
                end
                if W(3) > 0 % discr from D2, cont from D1
                    wtmp = Ddisc2.DetermY;
                    dtmp(discreten2) = Distribution();
                    for i=1:discreten2
                        dtmp(i) = pullBelowToX(Dcont1, Ddisc2.DetermX(i));
                    end
                    D(3) = wtmp * dtmp;                    
                end
                 % continous part: (both need a continuous part to let the following have a sense.)
                if W(4) > 0
                    ContPdfX = Dcont1.PdfX;
                    ContPdfY = Dcont1.PdfY .* Dcont2.CdfY + Dcont1.CdfY .* Dcont2.PdfY;                                        
                    D(4) = Distribution(DType.numeric, ContPdfX, ContPdfY, [],[], resultNonNormalized);
                end
                obj = W * D;
                
                % update support:
                obj.SupportLeft = sl;
                obj.SupportRight = sr;
                obj = obj.updateCdf;
                
                obj = obj.normalize();
            end
%             if ~done  
%                 
%                 discreten1 = length(D1.DetermX);
%                 discreten2 = length(D2.DetermX);
%                 discrprob1 = sum(D1.DetermY);
%                 discrprob2 = sum(D2.DetermY);
%                 W = [discrprob1 * discrprob2, discrprob1 * (1-discrprob2), (1-discrprob1)*discrprob2, (1-discrprob1)*(1-discrprob2)];
%                 D(4) = Distribution();
%                 % discrete part:
%                 if W(1) > 0
%                     newDetX = sort(unique([D1.DetermX, D2.DetermX]));
% %                     newDetY = newDetX * 0;
%                     DetX1    = [-Inf, D1.DetermX];
%                     cumDetY1 = [0,    cumsum(D1.DetermY)] / discrprob1;
%                     DetX2    = [-Inf, D2.DetermX];
%                     cumDetY2 = [0,    cumsum(D2.DetermY)] / discrprob2;
%                     cumMaxXY = newDetX * 0;
%                     for i=1:length(cumMaxXY)
%                         i1 = find(DetX1 <= newDetX(i), 1, 'last');
%                         i2 = find(DetX2 <= newDetX(i), 1, 'last');
%                         cumMaxXY(i) = cumDetY1(i1) * cumDetY2(i2);
%                     end
%                     newDetY = diff([0 cumMaxXY]);
%                     D(1) = Distribution(DType.deterministic, newDetX, newDetY);
%                 else
%                     D(1) = Distribution(DType.deterministic, 0,1);
%                 end
%                 % hybrid part:
%                 if W(2) > 0 % discr from D1, cont from D2
%                     wtmp = D1.DetermY / discrprob1;
%                     dtmp(discreten1) = Distribution();
%                     for i=1:discreten1
%                         dtmp(i) = pullBelowToX(D2, D1.DetermX(i));
%                     end
%                     D(2) = wtmp * dtmp;
%                 end
%                 if W(3) > 0 % discr from D2, cont from D1
%                     wtmp = D2.DetermY / discrprob2;
%                     dtmp(discreten2) = Distribution();
%                     for i=1:discreten2
%                         dtmp(i) = pullBelowToX(D1, D2.DetermX(i));
%                     end
%                     D(3) = wtmp * dtmp;                    
%                 end
%                  % continous part: (both need a continuous part to let the following have a sense.)
%                 if W(4) > 0
%                     ContPdfX = D1.PdfX;
% %                     ContPdfY = D1.PdfY .* D2.CdfY + D1.CdfY .* D2.PdfY;                    
%                     ContPdfY = D1.PdfY .* D2.CdfY / (1-discrprob1) + D1.CdfY .* D2.PdfY / (1-discrprob2);                    
%                     D(4) = Distribution(DType.numeric, ContPdfX, ContPdfY);
%                 end
%                 obj = W * D;
%                 
%                 % update support:
%                 obj.SupportLeft = sl;
%                 obj.SupportRight = sr;
%                 obj = obj.updateCdf;
%                 
%                 obj = obj.normalize();
%             end
        end
        %% --------------------------------------------------------------- %
        function obj = minDist(D1, D2)
            % Distribution of min(D1, D2) for independent distributions D1, D2
            if nargin < 2 % D1 hopefully contains a vector in this case
                if length(D1) == 1
                    obj = D1; 
                    return
                end                
                ihalf = ceil(length(D1) / 2);
                D2 = D1(ihalf+1:end);
                D1 = D1(1:ihalf);
            end
            
            if Distribution.DebugMode
                fprintf('DEBUG: minDist([%s], [%s])\n', sprintf('%d ', D1.ID), sprintf('%d ', D2.ID));
            end
            
            % TODO: elementwise?!
            
            if length(D1) > 1
                %recursively:
                D1 = minDist(D1);
            end
            if length(D2) > 1
                %recursively:
                D2 = minDist(D2);
            end
            
            obj = - maxDist(-D1, -D2);            
        end
        %% --------------------------------------------------------------- %
        function obj = inverseDist(obj)
            % determines the distribution of the inverse RV.
            % only possible for strictly positive RVs (in the continous domain).
            
            if Distribution.DebugMode
                fprintf('DEBUG: inverseDist([%s])\n', sprintf('%d ', obj.ID));
            end
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                res = Distribution.initDist(size(obj));
                for oi=1:numel(res)
                    res(oi) = obj(oi).inverseDist();
                end
                obj = res;
                return;
            end
            %
            
            % continous part:
            x = obj.PdfX;
            y = obj.PdfY;
            
            assert( all( y(x<=0) == 0 ) );
            
            newx = x;
            newy = obj.interpPDF( 1 ./ newx ) ./ (newx).^2;
            newy(newx<= 0) = 0;
            
            % discrete part
            discx = obj.DetermX;
            discy = obj.DetermY;
            
            assert(all(discx(discy>0) ~= 0));
            newdiscx = 1./discx;
            newdiscy = discy;
            
            if all(newy == 0)
                % pure discrete Distribution:
                obj = Distribution(DType.deterministic, newdiscx, newdiscy);
            else
                % general numeric Distribution:
                obj = Distribution(DType.numeric, newx, newy, newdiscx, newdiscy);
            end
        end
        %% --------------------------------------------------------------- %
        function samples = samples(obj, msamples, nsamples)
            % numerically draw nsamples samples from the Distribution:
%             [uy, ia, ~] = unique(obj.CdfY);
%             ux = obj.CdfX(ia);
            if Distribution.DebugMode
                fprintf('DEBUG: samples(%d)\n', obj.ID);
            end

            if obj.Type == DType.invalid 
                samples = NaN(msamples, nsamples); %return invalid
                return;
            end

            if obj.NonNormalized && Distribution.DistributionWarningLvl > 1
                warning('Creating Samples from nonNormalized distribution');
            end

            if nargin < 3
                nsamples = msamples;
                msamples = 1; % to get a row vector as before.
            end
            % Distribution specific sampling might be provided:
            switch (obj.Type)
%                 case DType.exp % REQUIRES STATISTICS TOOLBOX
%                     samples = random('exp', 1/obj.Params, msamples, nsamples);
                otherwise
                    ry = rand(1, msamples*nsamples);
                    % version 1: joint consideration:
%                     [ux, uy] = obj.ucdf();
%                     samples = interp1(uy, ux, ry);
                    
                    % version 2: split continuous and discrete part and treat separately:
                    [Ddisc, Dcont, probs] = getDiscreteAndContinuousPart(obj);
                    inddisc = ry <= probs(1);
                    indcont = ~inddisc;
                    
                    % discrete part
                    if any(inddisc)
                        s = zeros(1, sum(inddisc));
                        cumY = cumsum(Ddisc.DetermY);
                        rydisc = ry(inddisc) / probs(1);
                        for i=1:length(s)
                            si = find(rydisc(i) <= cumY , 1, 'first');
                            s(i) = Ddisc.DetermX(si);
                        end
                        samples(inddisc) = s;
                    end
                    
                    % continuous part:
                    if any(indcont)
                        [ux, uy] = Dcont.ucdf();
                        rycont = (ry(indcont)-probs(1)) / probs(2);
                        samples(indcont)= interp1(uy, ux, rycont );
                    end
                    
                    % reshaping:
                    samples = reshape(samples, [msamples, nsamples]);
                    
            end
        end
        %% --------------------------------------------------------------- %
        function val = percentile(obj, p)
            % returns the (100*p)'s percentile of a distribution.
            % I.e., p is given as a value between 0 and 1.
            if Distribution.DebugMode
                fprintf('DEBUG: percentile(%d, %g)\n', obj.ID, p);
            end
            assert(all(0 < p) && all(p < 1));
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                val = nan(numel(obj), numel(p));
                for oi=1:numel(obj)
                    val(oi, :) = obj(oi).percentile(p);
                end
                return;
            end
            %
            
            if obj.Type == DType.invalid
                val = nan; %return invalid
                return;
            end
            
            if obj.NonNormalized && Distribution.DistributionWarningLvl > 1
                warning('Calculating percentile from nonnormalized distribution');
            end
            
%             uy = obj.CdfY;
%             ux = obj.CdfX;
%             [uy, ia, ~] = unique(uy);
%             ux = ux(ia);    
            
            % instable version:
%             [ux, uy] = obj.ucdf();
%             val = interp1(uy, ux, p);
            
            % more stable version:
            val = nan(size(p));
            for pi = 1:length(p)
                pp = p(pi);
                ii = find(obj.CdfY > pp, 1, 'first');
                if isempty(ii) 
                    % no value is > pp:
                    val(pi) = +Inf;
%                     val(pi) = nan;
                elseif ii==1
                    % all values are > pp:
                    val(pi) = -Inf;
                else
                    % if obj.CdfX(ii-1) < obj.CdfX(ii) % given
                    val(pi) = interp1( obj.CdfY(ii-1:ii), obj.CdfX(ii-1:ii), pp);
                end
            end
        end
        %% --------------------------------------------------------------- %
        function val = median(obj)
            if Distribution.DebugMode
                fprintf('DEBUG: median([%s])\n', sprintf('%d ', obj.ID));
            end
            val = obj.percentile(0.5);
        end
        %% --------------------------------------------------------------- %
        function [ux, uy] = ucdf(obj) % unique cdf
            if Distribution.DebugMode
                fprintf('DEBUG: ucdf(%d)\n', obj.ID);
            end
            if obj.Type == DType.invalid
                ux = Distribution.DefaultPdfX; %return invalid
                uy = nan(size(ux));
                return;
            end
            
            uy = obj.CdfY;
            ux = obj.CdfX;
            [uy, ia, ~] = unique(uy);
            ia = sort(ia);              % don't mess up the original order
            ux = ux(ia);
            ia = ia';
            ind = diff([0 ia]) > 1;     % find each position with multiple entries
            fak = 2;
            if any(ind)
                ind = find(ind);                
                uy = [uy uy(ind-1)+fak*eps(uy(ind-1))]; % add some slightly differing values
                ux = [ux ux(ind)];
                if ia(end) < length(obj.CdfY)
                    uy = [uy obj.CdfY(end)+fak*eps(obj.CdfY(end))];
                    ux = [ux obj.CdfX(end)];
                end
                [uy, ui] = sort(uy);
                ux = ux(ui);
                %if this procedure introduced new duplicates, those can be safely removed:
                [uy, ia, ~] = unique(uy);
                ux = ux(ia);
            end            
        end
        %% --------------------------------------------------------------- %
        function valid = isvalid(obj)
%             if nargin == 
                valid = ~([obj.Type] == DType.invalid);
%             end
        end
        %% --------------------------------------------------------------- %
        function [h, ph, dh, ch] = plot(varargin)
%             if Distribution.DebugMode
%                 fprintf('DEBUG: plot()\n');
%             end
            ph = []; 
            dh = [];
            ch = [];
            h  = nan;
%             if nargin == 1
            if ~ishandle(varargin{1}(1))
                ax = gca;
                obj = varargin{1};
                vaii = 2;
            else                
                ax = varargin{1};
                obj = varargin{2};
                vaii = 3;
            end     
            
            if Distribution.DebugMode
                fprintf('DEBUG: plot([%s])\n', sprintf('%d ', obj.ID));
            end
            
            if numel(obj) > 1
                % recursive plotting of the whole vector
                hold on;
                h = nan(size(obj));
                ph = nan(size(obj));
                dh = nan(size(obj));
                ch = nan(size(obj));
                for i=1:numel(obj)
                    if nargin > 2
                        [hnew, phnew, dhnew, chnew] = plot(ax, obj(i), varargin{vaii:end});
                    else
                        [hnew, phnew, dhnew, chnew] = plot(ax, obj(i));
                    end
                    h (i) =  hnew;
                    if ~isempty(phnew)
                        ph(i) = phnew;
                    end
                    if ~isempty(dhnew)
                        dh(i) = dhnew;
                    end
                    if ~isempty(chnew)
                        ch(i) = chnew;
                    end
                end
                % we are already done here.
                return
            end           
            
            if obj.Type == DType.invalid
                return
            end
            
            colors = get(ax, 'ColorOrder');
            ci = get(ax, 'ColorOrderIndex');
            if ci > size(colors,1)
                ci = 1;
            end
%             if nargin == 1
% %                 obj = varargin{1};
%                 ph = plot(obj.PdfX,obj.PdfY);
%                 hold on;
%                 dh = plot(obj.DetermX, obj.DetermY, 'x', 'Color', colors(ci,:));
%                 ch = plot(obj.CdfX,obj.CdfY, '--', 'Color', colors(ci,:));
%             elseif nargin == 2
% %                 ax = varargin{1};
% %                 obj = varargin{2};
%                 ph = plot(ax,obj.PdfX,obj.PdfY);
%                 hold on;
%                 dh = plot(ax, obj.DetermX, obj.DetermY, 'x', 'Color', colors(ci,:));
%                 ch = plot(ax,obj.CdfX,obj.CdfY, '--', 'Color', colors(ci,:));
%             elseif nargin == 3
%                 warning('To use extra parameters use plot(ax, D, options)');
%             elseif nargin > 3
% %                 ax = varargin{1};
% %                 obj = varargin{2};
%                 ph = plot(ax,obj.PdfX,obj.PdfY,varargin{3:end});
%                 hold on;  
%                 dh = plot(ax, obj.DetermX, obj.DetermY, 'x',varargin{3:end});              
%                 ch = plot(ax,obj.CdfX,obj.CdfY, '--',varargin{3:end});
%             end
            
            % look for any Distribution specific options:
            opts = {};
            if nargin >= 3
                opts = varargin(3:end);
            end
            [plotpdf, opts] = parseoptions(opts, 'pdf', 'bool', 1);
            [plotcdf, opts] = parseoptions(opts, 'cdf', 'bool', 1);
            [ccdf, opts]    = parseoptions(opts, 'ccdf', 'bool', 0);
            [integratediscrete, opts] = parseoptions(opts, 'integratediscrete', 'bool', 0);
            [integratediscrete, opts] = parseoptions(opts, 'intdisc', 'bool', integratediscrete);
            defaultcdfline = '--';
            if ~plotpdf
                defaultcdfline = '-';
            end
            [cdfLine, opts] = parseoptions(opts, 'cdfline', 'char', defaultcdfline);
            
            [color, opts] = parseoptions(opts, 'Color', 'double', colors(ci,:));
            

            px = obj.PdfX;
            py = obj.PdfY;
%             integratediscrete=true;
            if integratediscrete
                for i = 1:length(obj.DetermX)
                    x = obj.DetermX(i);
                    y = obj.DetermY(i);
                    y0 = 0;
                    if x <= px(1)
                        px = [ x * ones(1,3), px ];
                        if x == px(1)
                            y0 = py(1); 
                        end
                        py = [ 0, y0+y, y0, py ];
                        continue
                    end
                    if x >= px(end)
                        px = [ px, x * ones(1,3) ];
                        if x == px(end)
                            y0 = py(end); 
                        end
                        py = [ py, y0, y0+y, 0];
                        continue
                    end
                    pos = find(px <= x, 1, 'last');
                    y0 = interp1( px(pos:pos+1), py(pos:pos+1), x );
                    px = [ px(1:pos), x * ones(1,3), px(pos+1:end) ];
                    py = [ py(1:pos), y0, y0+y, y0, py(pos+1:end) ];
                end
            end
            
            hold on;
            % if ccdf should be plotted:
            ycdfdata = obj.CdfY;
            if ccdf
                ycdfdata = 1-ycdfdata;
            end
            if isempty(opts)
                if plotpdf
                    ph = plot(ax, px, py, 'Color', color);
                    if ~integratediscrete
%                         dh = plot(ax, obj.DetermX, obj.DetermY, 'x', 'Color', colors(ci,:));
                        dh = stem(ax, obj.DetermX, obj.DetermY, 'Color', color);
                    end
                end                
                if plotcdf                    
                    ch = plot(ax,obj.CdfX,ycdfdata, cdfLine, 'Color', color);
                end
            else            
%                 warning('Experimental');
                
                if plotpdf
                    ph = plot(ax, px, py, opts{:}, 'Color', color);
                    if ~integratediscrete
                        dh = stem(ax, obj.DetermX, obj.DetermY, opts{:}, 'Color', color);              
                    end
                end
                if plotcdf
%                     ch = plot(ax,obj.CdfX,obj.CdfY, cdfLine, opts{:}, 'Color', color);
                    ch = plot(ax,obj.CdfX,ycdfdata, opts{:}, 'Color', color);
                end
            end
            
%             if isempty(ph)
%                 ph = dh;
%             end
%             if isempty(dh)
%                 dh = ph;
%             end
            
            % return representive handle:
            if ~plotpdf
                h = ch;
            else
                if obj.Type == DType.deterministic
                    h = dh;
                else
                    h = ph;                    
                end
            end
            
            % remove other plots from legend:
            if ~isempty(ph) && ph~=h
                set(ph, 'HandleVisibility', 'off');
            end
            if ~isempty(ch) && ch~=h
                set(ch, 'HandleVisibility', 'off');
            end
            if ~isempty(dh) && dh~=h
                set(dh, 'HandleVisibility', 'off');
            end
            
            % don't skip colors
            set(ax, 'ColorOrderIndex', ci+1);
%             h = [ph; dh; ch];
%             xlabel('Time')
%             ylabel('pdf')
        end
        %% --------------------------------------------------------------- %
        function print(obj)
            if length(obj) > 1
                fprintf('Array of Distributions:\n');
                for i=1:length(obj)
                    fprintf('Element %d:\n', i);
                    obj(i).print();
                end
                return
            end
            
            if Distribution.DebugMode
                fprintf('DEBUG: print(%d)\n', obj.ID);
            end
            
%             fprintf('Distribution (%s)\n', char(obj.Type));
            fprintf('Distribution %s\n', obj.shortname());
            fprintf('ID: %4d    Created: %s\n', obj.ID, datestr(obj.created));
            fprintf('NonNormalized: %d\n', obj.NonNormalized);
            if ~isempty(obj.Params)
                fprintf('Parameter: \n');
                if length(obj.ParamNames) == length(obj.Params)
                    for i=1:length(obj.Params)
                        fprintf('%g\t%s\n', obj.Params(i), obj.ParamNames{i});
                    end
                else
                    fprintf('(Parameter - Name - mismatch)\n');
                    fprintf('%g\t', obj.Params);
                end
                fprintf('\n');
            end
            if ~isempty(obj.DetermX)
                fprintf('Deterministic Part:\n');
                nl = 10;
                nd = length(obj.DetermX);
                for i=1:ceil(nd / nl)
                    ind = 1+(i-1)*nl:min(nd, i*nl);
                    fprintf('%9g\t', obj.DetermX(ind));
                    fprintf('\n');
                    fprintf('%9g\t', obj.DetermY(ind));
                    fprintf('\n\n');
                end
            end
        end
        %% --------------------------------------------------------------- %
        function str = shortname(obj)
            
            % Handling Arrays:
            if numel(obj) > 1
                str = cell(size(obj));
                for oi=1:numel(obj)
                    str{oi} = obj(oi).shortname();
                end
                return
            end
            %
            
            switch obj.Type
                case DType.exp
                    str = sprintf('Exp(%g)', obj.Params(1));
                case DType.uniform
                    str = sprintf('Uniform([%g, %g])', obj.Params(1:2));
                case DType.deterministic                    
                    valstr = sprintf('%g ', obj.DetermX);
                    if length(valstr) > 1
                        valstr(end) = []; %remove last space
                    end
                    if length(obj.DetermX) == 1
                        str = sprintf('Deterministic(%s)', valstr);
                    else
                        str = sprintf('Discrete ([%s])', valstr); 
                    end
                    
                otherwise
                    str = char(obj.Type);                    
            end            
        end
        %% --------------------------------------------------------------- %
        function dx = getdx(obj)
            if Distribution.DebugMode
                fprintf('DEBUG: getdx(%d)\n', obj.ID);
            end
            dx = diff(obj.PdfX(1:2));
        end
        %% --------------------------------------------------------------- %
        function axlim = getlimits(obj, fac, pdfcdf)
            if Distribution.DebugMode
                fprintf('DEBUG: getlimits([%s])\n', sprintf('%d ', obj.ID));
            end
            if nargin < 2 || isempty(fac)
                fac = 0;
            end
            if nargin < 3 || isempty(pdfcdf)
                pdfcdf = 'pdf';
            end
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                axlim = nan(numel(obj), 4);
                for oi=1:numel(obj)
                    axlim(oi,:) = obj(oi).getlimits(fac, pdfcdf);
                end
                return;
            end
            %
            
            left  = obj.percentile(0.01);
            right = obj.percentile(0.99);
            if strcmpi(pdfcdf, 'pdf')
                top   = max( [obj.PdfY, obj.DetermY] );
            else
                top   = 1;
            end
            bottom = 0;
            w = right - left;
            h = top - bottom;
            axlim = [ left - fac*w, right + fac*w, bottom-fac*h, top+fac*h ];
        end
        %% --------------------------------------------------------------- %
        function Y = interpPDF(obj, X)
            if Distribution.DebugMode
                fprintf('DEBUG: interpPDF([%s], x)\n', sprintf('%d ', obj.ID));
            end
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                Y = nan(numel(obj), numel(X));
                for oi=1:numel(obj)
                    Y(oi,:) = obj(oi).interpPDF(X);
                end
                return;
            end
            %
            
            if obj.Type == DType.invalid 
                Y = NaN(size(X)); %return invalid
                return;
            end
            Y = zeros(size(X));
            valid = X >= obj.PdfX(1) & X <= obj.PdfX(end);
            Y(valid) = interp1(obj.PdfX, obj.PdfY, X(valid));
        end
        %% --------------------------------------------------------------- %
        function Y = interpCDF(obj, X)
            if Distribution.DebugMode
                fprintf('DEBUG: interpCDF([%s], x)\n', sprintf('%d ', obj.ID));
            end
            
            % Handling of arrays of Distributions:
            if numel(obj) > 1
                Y = nan(numel(obj), numel(X));
                for oi=1:numel(obj)
                    Y(oi,:) = obj(oi).interpCDF(X);
                end
                return;
            end
            %
            
            if obj.Type == DType.invalid 
                Y = NaN(size(X)); %return invalid
                return;
            end
            Y = zeros(size(X));
            valid = X >= obj.CdfX(1) & X <= obj.CdfX(end);
%             [ux, uy] = ucdf(obj); 
            Y(valid) = interp1(obj.CdfX, obj.CdfY, X(valid));
            
            if obj.NonNormalized
                Y(X>obj.CdfX(end)) = max(obj.CdfX);                
            else
                Y(X>obj.CdfX(end)) = 1;
            end
        end
        % --------------------------------------------------------------- %
        function splitD = split(obj, p)
            % splits a process into multiple processes according to the probabilities given in p
            if nargin < 2 || isempty(p)
                if Distribution.DistributionWarningLvl > 1
                    warning('no p specified in split');
                end
                p = [0.5 0.5];
            end
            if Distribution.DebugMode
                fprintf('DEBUG: split(%d, [%s])\n', obj.ID, sprintf('%g ', p));
            end
            nsplit = length(p);
            
            % NO ARRAY IMPLEMENTATION YET
            
            nmax = 30;            
            if length(obj) == nmax
                % assume that Dsum is precalculated:
                Dsum = obj;
            else
                Dsum = getSumDists(obj, nmax);
            end
            exponent = (1:nmax)-1;
            splitD(nsplit) = Distribution();
            for si = 1:nsplit
                w = (1-p(si)) .^ exponent * p(si);
                splitD(si) = w * Dsum;
            end            
        end
        %% --------------------------------------------------------------- %
        function obj = pullNegativeToZero(obj)
            if Distribution.DebugMode
                fprintf('DEBUG: pullNegativeToZero(%d)\n', obj.ID);
            end
             % NO ARRAY IMPLEMENTATION YET
            
%             obsolete, retained as a special case of the more general new function:
            obj = pullBelowToX(obj, 0);
%             ineg = find(obj.CdfX < 0, 1, 'last');
%             % probability of being negative:
%             prob = obj.CdfY(ineg); % should include discrete and continious parts
%             if isempty(prob)
%                 prob = 0;
%             end
%             % remove negative part:
%             obj.PdfY(1:ineg) = 0;
%             obj.CdfY(1:ineg) = 0;
%             obj.DetermY(obj.DetermX < 0 ) = 0;
%             % add removed part to a deterministic zero.
%             izero = find(obj.DetermX == 0 );
%             if ~isempty(izero)
%                 obj.DetermY(izero(1)) = obj.DetermY(izero(1))+prob;
%             else
%                 obj.DetermX = [obj.DetermX 0];
%                 obj.DetermY = [obj.DetermY prob];
%             end
        end
        %% --------------------------------------------------------------- %
        function obj = pullBelowToX(obj, X)
             % NO ARRAY IMPLEMENTATION YET
             
            if nargin < 2
                X = 0;
            end
            
            if Distribution.DebugMode
                fprintf('DEBUG: pullBelowToX(%d,%g)\n', obj.ID, X);
            end
            
            if obj.Type == DType.invalid 
                obj = Distribution(); %return invalid
                return;
            end
            ibelow = find(obj.CdfX < X, 1, 'last');
            % probability of being negative:
            prob = obj.CdfY(ibelow); % should include discrete and continious parts
            if isempty(prob)
                prob = 0;
            end
            % remove negative part:
            % continuous:
            obj.PdfY(1:ibelow) = 0;
            obj.CdfY(1:ibelow) = 0;
            % discrete:
%             obj.DetermY(obj.DetermX < X ) = 0;
            obj.DetermY(obj.DetermX < X ) = [];
            obj.DetermX(obj.DetermX < X ) = [];
            % add removed part to a deterministic X.
            iX = find(obj.DetermX == X );
            if ~isempty(iX)
                obj.DetermY(iX(1)) = obj.DetermY(iX(1))+prob;
            else
                obj.DetermX = [obj.DetermX X];
                obj.DetermY = [obj.DetermY prob];
            end
        end
        %% --------------------------------------------------------------- %
        function [Ddisc, Dcont, probs] = getDiscreteAndContinuousPart(obj)
            if Distribution.DebugMode
                fprintf('DEBUG: getDiscreteAndContinuousPart(%d)\n', obj.ID);
            end
            if obj.Type == DType.invalid 
                Ddisc = Distribution(); %return invalid
                Dcont = Distribution();
                probs = [NaN NaN];
                return;
            end
            
            if obj.NonNormalized && Distribution.DistributionWarningLvl > 1
                warning('getDiscreteAndContinuousPart from nonnormalized Distribution!');
            end
            discPartProb = sum(obj.DetermY);
            if obj.NonNormalized
                probs = [discPartProb, obj.CdfY(end)-discPartProb];
            else
                probs = [discPartProb, 1-discPartProb];
            end
            
            if probs(1) == 0 %discPartProb == 0
                Dcont = obj;
                Ddisc = Distribution(DType.invalid);
            elseif probs(2) == 0 % discPartProb == 1
                Dcont = Distribution(DType.invalid);
                Ddisc = obj;
            else
                Dcont = Distribution(DType.numeric, obj.PdfX, obj.PdfY / probs(2), [], [], obj.NonNormalized );
                Ddisc = Distribution(DType.deterministic, obj.DetermX, obj.DetermY / probs(1), obj.NonNormalized );
            end
        end
        %% --------------------------------------------------------------- %
        function [DEqualTo, DRemaining, probs] = ExtractDiscEqualTo(obj, EqualTo)
            if Distribution.DebugMode
                fprintf('DEBUG: ExtractDiscEqualTo(%d)\n', obj.ID);
            end
            % invalid:
            if obj.Type == DType.invalid
                DEqualTo = Distribution(); %return invalid
                DRemaining = Distribution();
                probs = [NaN NaN];
                return;
            end
            if obj.NonNormalized && Distribution.DistributionWarningLvl > 1
                warning('ExtractDiscEqualTo from nonnormalized Distribution!');
            end
            DEqualTo = Distribution(DType.deterministic, EqualTo, 1);
            DRemaining = obj;
            
            indEqualTo = find(obj.DetermX == EqualTo);
            if ~isempty(indEqualTo) && sum(obj.DetermY(indEqualTo)) > 0
                p = sum(obj.DetermY(indEqualTo));
                probs=[p, 1-p];
                
                if p == 1
                    % nothing left. Return invalid
                    DRemaining = Distribution();
                    return;
                end
                
                % remove the respective part:
                DRemaining.DetermX(indEqualTo) = [];
                DRemaining.DetermY(indEqualTo) = [];
                
                % rescale the remaining
                DRemaining.DetermY = DRemaining.DetermY / (1-p);
                DRemaining.PdfY = DRemaining.PdfY / (1-p);
                % either update cdf:
                DRemaining = DRemaining.updateCdf();
                % or create a new distribution?!
                if Distribution.DistributionWarningLvl > 1
                    warning('rethink this... (ExtractDiscEqualTo(...))');
                end
            else
                probs=[0 1];
            end            
        end        
        %% --------------------------------------------------------------- %
        function [DNeg, DZero, DPos, probs] = splitNegZeroPos(obj, thresh)
            if nargin < 2 || isempty(thresh)
                thresh = 0;
            end
            probs = [0 0 0];
            [DZero, DNonZero, probszerononzero] = ExtractDiscEqualTo(obj, thresh);
            probs(2) = probszerononzero(1);

            if DNonZero.isvalid
                %
                probs(1) = DNonZero.interpCDF(thresh);
                probs(3) = 1-probs(1);
                %
                discindp = DNonZero.DetermX > thresh;
                contindp = DNonZero.PdfX > thresh;
                discindn = DNonZero.DetermX < thresh;
                contindn = DNonZero.PdfX < thresh;
                DNeg = DNonZero;
                DPos = DNonZero;
                
                if probs(1) == 0
                    DNeg = Distribution();
                else
                    DNeg.DetermX(~discindn) = [];
                    DNeg.DetermY(~discindn) = [];
                    DNeg.PdfY(~contindn)    = 0;
                    DNeg.DetermY = DNeg.DetermY / probs(1);
                    DNeg.PdfY    = DNeg.PdfY    / probs(1);
                    DNeg = DNeg.updateCdf();
                end

                if probs(3) == 0
                    DPos = Distribution();
                else
                    DPos.DetermX(~discindp) = [];
                    DPos.DetermY(~discindp) = [];
                    DPos.PdfY(~contindp)    = 0;
                    DPos.DetermY = DPos.DetermY / probs(3);
                    DPos.PdfY    = DPos.PdfY    / probs(3);
                    DPos = DPos.updateCdf();
                end

                %               
                probs([1,3]) = probs([1,3]) * probszerononzero(2);
            else
                DNeg = Distribution();
                DPos = Distribution();
            end
        end
        %% --------------------------------------------------------------- %
        function [DNonNeg, DNeg, probs] = getNonNegativePart(obj)
            [DNeg, DZero, DPos, probs1] = splitNegZeroPos(obj);
            probs = [1-probs1(1), probs1(1)];
            if probs1(1) < 1
                DNonNeg = probs1(2:3) / sum(probs1(2:3)) * [DZero, DPos];
            else
                DNonNeg = Distribution();
            end
        end
        %% --------------------------------------------------------------- %
        function [DPos, DNonPos, probs] = getPositivePart(obj)
            [DNeg, DZero, DPos, probs1] = splitNegZeroPos(obj);
            probs = [probs1(3), 1-probs1(3)];
            if probs1(1) < 1
                DNonPos = probs1(1:2) / sum(probs1(1:2)) * [DNeg, DZero];
            else
                DNonPos = Distribution();
            end
        end
        %% --------------------------------------------------------------- %
        function d = distance(D1, D2, whichnorm, pdfcdf, domain)
            if nargin < 3 || isempty(whichnorm)
                whichnorm = 2;
            end
            if nargin < 4 || isempty(pdfcdf)
                pdfcdf = 'pdf';
            end
            if nargin < 5 || isempty(domain)
                domain = [-inf, +inf];
            end
            
            if Distribution.DebugMode
                fprintf('DEBUG: distance(%d, %d, %g, %s)\n', D1.ID, D2.ID, whichnorm, pdfcdf);
            end
            
            if D1.Type == DType.invalid || D2.Type == DType.invalid
                d = Inf;
                return
            end
            
            %only implemented for default domains.
            assert(D1.isdefault());
            assert(D2.isdefault());            
            
            % take care of deterministic part:
            nx1 = length(D1.DetermX);
            [commonDetX, ~,ic] = unique( [D1.DetermX, D2.DetermX] );
            commonDetXY1 = zeros(size(commonDetX));
            commonDetXY2 = zeros(size(commonDetX));
            
            commonDetXY1(ic(1:nx1)) = D1.DetermY;
            commonDetXY2(ic((nx1+1):end)) = D2.DetermY;
            
            if strcmpi(pdfcdf, 'pdf')
                ContY1 = D1.PdfY;
                ContY2 = D2.PdfY;
                Disc1 = commonDetXY1;
                Disc2 = commonDetXY2;
            else
                ContY1 = D1.CdfY;
                ContY2 = D2.CdfY;
                Disc1 = commonDetXY1 * 0; % ignore this in this case
                Disc2 = commonDetXY2 * 0; % ignore this in this case
            end
            
            % apply the domain:
            ind = D1.PdfX < domain(1) | D1.PdfX > domain(2);
            ContY1( ind ) = [];
            ContY2( ind ) = [];
            
            ind = commonDetX < domain(1) | commonDetX > domain(2);
            Disc1( ind ) = [];
            Disc2( ind ) = [];
            
            switch whichnorm
                case 2
                    d = sum((ContY1-ContY2).^2);
                    d = d + sum( (Disc1 - Disc2).^2 );
                    d = sqrt(d);
                otherwise
                    if Distribution.EnableDistributionWarnings 
                        warning('unknown norm');
                    end
                    d=nan;
            end
        end
        % --------------------------------------------------------------- %
        
        % --------------------------------------------------------------- %
        
        % --------------------------------------------------------------- %
        
        % --------------------------------------------------------------- %
    end
    
    % ------------------------------------------------------------------- %
    % STATIC METHODS ---------------------------------------------------- %
    % ------------------------------------------------------------------- %
     methods (Static)
%         function SetDefaultPdfX(NewPdfX)
%             Distribution.DefaultPdfX = NewPdfX;
%         end
        function [getdata, getDomainChanged] = DefaultPdfX(setdata)
            persistent persDefaultPdfX;
            persistent DomainChanged;
            defaultdomain = linspace(-1000,2000,60001);
            if nargin % can be changed:
                if ischar(setdata) && strcmpi(setdata, 'default')
                    persDefaultPdfX = defaultdomain; 
                    DomainChanged = false;
                else
                    persDefaultPdfX = setdata; 
                    DomainChanged = true;
                end
            end
            if isempty(persDefaultPdfX) % default value if not set at the begining
                persDefaultPdfX = defaultdomain; 
                DomainChanged = false;
            end
            getdata = persDefaultPdfX; 
            getDomainChanged = DomainChanged;
        end
        %% -------------------------------------------- %
        function out = currentID()
            persistent id;
            if isempty(id)
                id = 1; 
            end
            out = id;
            id = id+1;
        end
        %% -------------------------------------------- %
        function dx = getdefaultdx()
            dx = diff( Distribution.DefaultPdfX(1:2) );
        end
        %% -------------------------------------------- %
        function D = exampleDist(name, randomvalues)
            
            if nargin < 2 || isempty(randomvalues)
                randomvalues = false;
            end
                
            switch lower(name)
                case 'discrete'
                    if randomvalues
                        N = 1 + ceil( rand() * 9 ); % 2..10 values
                        discx = rand(1,N) * 20;     % values between 0..20
                        discy = rand(1,N);
                        discy = discy / sum(discy);
                    else
                        discx = [0 3 7 9 11 12];
                        discy = [10 12 7 4 2 1] / 36;
                    end
                    D = Distribution(DType.deterministic, discx, discy);
                case 'continuous'
                    % mixed gaussian
                    if randomvalues
                        N = 1 + ceil( rand() * 3 ); % 2..4 values
                        mus = rand(1,N) * 20;       % values between 0..20
                        sigmas = 1 + rand(1,N) * 9; % values between 1..10
                        ps  = rand(1,N);
                        ps  = ps / sum(ps);
                    else
                        N = 2;
                        mus = [2 7];
                        sigmas = [3 5];
                        ps = [0.4 0.6];
                    end
                    Dtmp(N) = Distribution();
                    for i=1:N
                        Dtmp(i) = Distribution(DType.gauss, mus(i), sigmas(i));
                    end
                    D = ps * Dtmp;
                case 'mixed'
                    if randomvalues
                        ps = rand();
                        ps = [ps, 1-ps];
                    else
                        ps = [0.5 0.5];
                    end
                    D = ps * ...
                        [Distribution.exampleDist('discrete', randomvalues), ...
                        Distribution.exampleDist('continuous', randomvalues)];
                otherwise
                    if Distribution.DistributionWarningLvl > 1
                        warning('unknown name "%s" for exampleDist()', name);
                    end
                    D = Distribution.exampleDist('mixed', randomvalues);
            end
        end
        %% --------------------------------------------------------------- %
        function obj = initDist(varargin)
            % initializes an array of invalid Distributions
            % either:
            %   initDist(sz) with sz being an array of dimensions
            % or
            %   initDist(dim1,dim2, ..., dimN)
            %
            % if a Distribution is added as first argument like 
            % either:
            %   initDist(D, sz) with sz being an array of dimensions
            % or
            %   initDist(D, dim1,dim2, ..., dimN)
            % then this will be the default Distribution for the
            % initialization
            
            if isa(varargin{1}, 'Distribution')
                D0 = varargin{1};
                offset = 1;
            else
                D0 = Distribution();
                offset = 0;
            end
            
            if nargin < 1+offset
                % return empty array of Distribution
                warning ('empty Distribution initialized');
                obj = Distribution();                
                obj(1) = [];
                return;
            elseif nargin == 1+offset
                sz   = varargin{1+offset};
                ndim = length(sz);
            else
                ndim = nargin - offset;
                sz   = ones(1, ndim);              
                for i=1:ndim
                    sz(i) = varargin{i+offset};
                end
            end
            
            obj = repmat(D0, sz(1), 1);
            % add the other dimensions:
            for i=2:ndim
                obj = repmat(obj, [ones(1, i-1), sz(i)] );
            end            
        end
        %% -------------------------------------------------------------------- %
        function D = getDistFromCdf(xcdf, ycdf)
            % numerical differentiation to get the pdf
            
            xpdf = Distribution.DefaultPdfX;
            yyy  = xpdf * 0;
            dx   = diff(xpdf(1:2));
            if length(xpdf) ~= length(xcdf) || any( xcdf ~= xpdf )
                below = xpdf < xcdf(1);
                above = xpdf > xcdf(end);
                valid = ~below & ~above;
                yyy(valid) = interp1(xcdf, ycdf, xpdf(valid));
%                 yyy(below) = ycdf(1);
%                 yyy(above) = ycdf(end);
                yyy(above) = 1;
            else % same domain
                yyy = ycdf;
            end            
            ypdf = [0, (yyy(2:end) - yyy(1:end-1))] / dx;
            
            D = Distribution(DType.numeric, xpdf, ypdf);
        end
    end
end