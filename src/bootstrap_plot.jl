cd(“$(homedir())/Desktop”)

srand(1234); using PyPlot;

Ts = [100:500:20100];
boot_mean = zeros(1,length(Ts)); boot_sdest = zeros(1,length(Ts));
boot_sd = zeros(1,length(Ts));
samp_mean = zeros(1,length(Ts)); samp_sd = zeros(1,length(Ts));

for t = 1:length(Ts) #For each sample size
    T = Ts[t];
    y = randn(T,1); #mean zero variance 5/3, supposed to be t distribution;
    samp_mean[t] = mean(y); samp_sd[t] = sqrt(var(y));
    sampmean = zeros(200,1);sampsd = zeros(200,1);
    for r = 1:200 #bootstrap replication
        index = floor(rand(1,T)*(T))+1;
        boot_y = y[index];
        sampmean[r] = mean(boot_y);
        sampsd[r] = sqrt(var(boot_y));
    end
    boot_mean[t] = mean(sampmean);
    boot_sdest[t] = mean(sampsd);
    # Bootstrap SD of sample mean;
    boot_sd[t] = sqrt(1/(200-1)*(sum((sampmean-boot_mean[t]).^2)));
end

grid(“on”)
#Sample
h = figure(1);
subplot(3,2,1)
plot(Ts, samp_mean’);
xlim(Ts[1],Ts[end]);
xlabel('sample size T');
ylabel('Mean of Ssample’);
title('Sample Mean With Different Sample Size of t_5 distribution');

subplot(3,2,2)
plot(Ts, samp_sd’);

xlim(Ts(1) Ts(end)]);
set(q, 'linewidth', 1, 'color', 'r', 'linestyle', '--');
xlabel('sample size T','fontsize',14);
ylabel('Standard Deviation of Bootstrap resample','fontsize',14);
title('Sample S.D. With Different Sample Size of t_5 distribution','fontsize',14);


%Bootstrap

subplot(3,2,3)
plot(Ts, boot_mean);
q = line([Ts(1) Ts(end)],[0 0]);
xlim([Ts(1) Ts(end)]);
set(q, 'linewidth', 1, 'color', 'r', 'linestyle', '--');
xlabel('sample size T','fontsize',14);
ylabel('Mean of Bootstrap resample','fontsize',14);
title('Bootstrap Resample Mean With Different Sample Size of t_5 distribution','fontsize',14);

subplot(3,2,4)
plot(Ts, boot_sdest);
q = line([Ts(1) Ts(end)],[sqrt(5/3) sqrt(5/3)]);
xlim([Ts(1) Ts(end)]);
set(q, 'linewidth', 1, 'color', 'r', 'linestyle', '--');
xlabel('sample size T','fontsize',14);
ylabel('Standard Deviation of Bootstrap resample','fontsize',14);
title('Bootstrap Resample S.D. With Different Sample Size of t_5 distribution','fontsize',14);

%deviation from sample to (resample from bootstrap);
subplot(3,2,5)
plot(Ts, abs(boot_mean - samp_mean));
q = line([Ts(1) Ts(end)],[0 0]);
xlim([Ts(1) Ts(end)]);
set(q, 'linewidth', 1, 'color', 'r', 'linestyle', '--');
xlabel('sample size T','fontsize',14);
ylabel('Difference','fontsize',14);
title('Absolute Deviation of Original Sample Mean from The Bootstrap Sample Mean','fontsize',14);

subplot(3,2,6)
plot(Ts, abs(boot_sdest - samp_sd));
q = line([Ts(1) Ts(end)],[0 0]);
xlim([Ts(1) Ts(end)]);
set(q, 'linewidth', 1, 'color', 'r', 'linestyle', '--');
xlabel('sample size T','fontsize',14);
ylabel('Difference','fontsize',14);
title('Absolute Deviation of Original Sample S.D. from The Bootstrap Sample S.D. ','fontsize',14);

savefig(“bootstrappen.pdf”);