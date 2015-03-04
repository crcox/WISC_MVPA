function summarize_jobs()
	jdat = loadjson('master.json');
    nlam = length(jdat.lambda);
    for i = 1:length(jdat.config)
        jobdir = sprintf('%03d',i-1);
        load(fullfile(jobdir,'fitObj.mat'));
        [nsubj, ncv] = size(fitObj);
        dps = zeros(nsubj, nlam);
        for ss = 1:nsubj
            tmp = zero(ncv, nlam);
            for cc = 1:ncv
                tmp(cc,:) = dprimeCV(fitObj.Y, fitObj.Yh, fitObj.testset);
            end
            dps(ss,:) = mean(tmp);
        end
        dp = mean(dps)
    end
end
