% testdistributionsystematic.m
%=================================================
% tests some operations for the distribution class
%=================================================
% domain = [-5, 20];
domain = [0.1, 15];
nsamples = 1e6;
% nonnegative = true;
supportmode = 'positive'; % nonnegative or anything else

types = {'discrete', 'continuous', 'mixed'};
nD = length(types);

clear sampleD;

sampleD(nD) = Distribution();
D(nD) = Distribution();
samples = zeros(nD, nsamples);
samples2 = zeros(nD, nsamples);
for i=1:nD
    fprintf('Create %s Distribution and samples... \n', types{i});
    D(i) = Distribution.exampleDist(types{i});
    switch supportmode
        case 'nonnegative'
            D(i) = D(i).getNonNegativePart();
        case 'positive'
            D(i) = D(i).getPositivePart();
    end
    samples(i,:) = D(i).samples(nsamples);    
    samples2(i,:) = D(i).samples(nsamples);    
    sampleD(i) = getDistFromSamples2( samples(i,:) );    
end

figure;
hold on;
sh = [nan nan];
for s=1:2
    sh(s) = subplot(2,1,s);
    if s==1
        notpdfcdf = 'cdf';
    else
        notpdfcdf = 'pdf';
    end
    for i=1:nD
        h = plot(gca, D(i), notpdfcdf, 'off', 'linestyle', '-');
        plot(gca, sampleD(i), notpdfcdf, 'off', 'linestyle', ':', 'color', get(h, 'color'));
    end
    xlim(domain);
end
lstr = repmat(types, 2,1); lstr = lstr(:).';
legend(lstr);
linkaxes(sh,'x');
subplot(2,1,1);
title('Example Functions');

%%
% remaining_real = @(x) rand(size(x)) .* x;
remaining_real = @(x) x - rand(size(x)) .* x;
remaining_dist = @(X) getRemainingTime(X);

inverse_real = @(x) 1./x;
inverse_dist = @(X) inverseDist(X);


fhandles = {@plus, @plus; ...}; % function pair real (must work elmentwise on vectors), Dist
    remaining_real, remaining_dist;...
    inverse_real, inverse_dist };
fhandles = { inverse_real, inverse_dist };

nf = size(fhandles, 1);

sbdim = ceil(sqrt( nD * (nD+1) ) );
sbdim = [ceil((nD * (nD+1)) / sbdim), sbdim ];

LineStyleSim = '--';

tic
for fi = 1:nf
    f = fhandles{fi, 1};
    fD = fhandles{fi, 2};
    
    fprintf('Testing %s ...\n', func2str(f));
       
    
    narg = nargin(f);
    
    figure;
    hold on;
    if narg == 1
        newdomain = sort(f(domain));
                
        %arg 1:
        clear funcD sampleD;
        funcD(nD) = Distribution();
        funcsampleD(nD) = Distribution();
        for i=1:nD
            fprintf('testing f(D%d)...\n', i);
            funcD(i) = fD( D(i) );
            funcSamples = f( samples(i,:) );
            funcsampleD(i) = getDistFromSamples2( funcSamples );
            
            subplot(nD, 2, 2*i-1);            
            title(sprintf('%s, %s', func2str(f), types{i}));
            plot(gca, funcD(i), 'cdf', 'off');
            plot(gca, funcsampleD(i), 'cdf', 'off', 'LineStyle', LineStyleSim);
            xlim(newdomain);
            subplot(nD, 2, 2*i);            
            plot(gca, funcD(i), 'pdf', 'off');
            plot(gca, funcsampleD(i), 'pdf', 'off', 'LineStyle', LineStyleSim);
            xlim(newdomain);
            drawnow;
        end
        
    elseif narg==2
        newdomain = sort(f(domain, domain));
        
        clear funcD sampleD;
        funcD(nD, nD) = Distribution();
        funcsampleD(nD, nD) = Distribution();
        
        
        ii = 1;
        for i=1:nD
            for j=i:nD
                fprintf('testing f(D%d, D%d)...\n', i, j);
                funcD(i,j) = fD( D(i), D(j) );
                funcSamples = f( samples(i,:), samples2(j,:) ); % independent samples
                funcsampleD(i,j) = getDistFromSamples2( funcSamples );

                subplot(sbdim(1), sbdim(2), ii); ii = ii+1;
                title(sprintf('%s, %s, %s', func2str(f), types{i}, types{j}));
                plot(gca, funcD(i,j), 'cdf', 'off');
                plot(gca, funcsampleD(i,j), 'cdf', 'off', 'LineStyle', LineStyleSim);
                xlim(newdomain);
                subplot(sbdim(1), sbdim(2), ii); ii = ii+1;    
                plot(gca, funcD(i,j), 'pdf', 'off');
                plot(gca, funcsampleD(i,j), 'pdf', 'off', 'LineStyle', LineStyleSim);
                xlim(newdomain);
                drawnow;
            end
        end
    else
        warning('strange number of arguments %d', narg);
    end
    legend('model', 'sim', 'location', 'southeast');

end
toc
fprintf('done\n');