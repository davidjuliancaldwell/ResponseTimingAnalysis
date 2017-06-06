level = 5;
wname = 'sym4';
npc = 'kais';
x = squeeze(epochedCortEco(:,:,10));
[x_sim, qual, npc] = wmspca(x ,level, wname, npc);

kp = 0;
chansInt = [2 10 40 42];

npc(1:4) = 0;

[x_sim, qual, npc] = wmspca(x ,level, wname, npc);

figure
for i = chansInt
    subplot(4,2,kp+1), plot(x (:,i)); axis tight;
    title(['Original signal ',num2str(i)])
    subplot(4,2,kp+2), plot(x_sim(:,i)); axis tight;
    title(['Simplified signal ',num2str(i)])
    kp = kp + 2;
end

%%

x = squeeze(epochedCortEco(:,2,10));
x_2 = squeeze(processedSig(:,2,10));

wt = modwt(x,10);
wtrec = zeros(size(wt));
wtrec(10:11,:) = wt(10:11,:);
y = imodwt(wtrec,'sym4');

figure
plot(x)
hold on
plot(y)
ylim([-3e-4 3e-4])

wt = modwt(x_2,10);
wtrec = zeros(size(wt));
wtrec(9:11,:) = wt(9:11,:);
y = imodwt(wtrec,'sym4');

figure
plot(x_2)
hold on
plot(y)
ylim([-3e-4 3e-4])

%%

x = squeeze(epochedCortEco(:,2,10));
x_2 = squeeze(processedSig(:,2,10));

wt = modwt(x,15);
wtrec(10:15,:) = wt(10:15,:);
y = imodwt(wtrec,'sym4');

figure
plot(x)
hold on
plot(y)
ylim([-3e-4 3e-4])

wt = modwt(x_2,15);
wtrec = zeros(size(wt));
wtrec(10:15,:) = wt(10:15,:);
y = imodwt(wtrec,'sym4');

figure
plot(x_2)
hold on
plot(y)
ylim([-3e-4 3e-4])
%%
x = squeeze(epochedCortEco(:,2,10));

wt = modwt(x,15);
wtrec = zeros(size(wt));
wtrec(9:15,:) = wt(9:15,:);
y = imodwt(wtrec,'sym4');

figure
plot(x)
hold on
plot(y)
ylim([-3e-4 3e-4])

figure
[f,P1] = spectralAnalysisComp(eco_fs,x);
plot(f,P1)

set(gca,'xscale','log');
set(gca,'yscale','log');

figure
[f,P1] = spectralAnalysisComp(eco_fs,y);
plot(f,P1)

set(gca,'xscale','log');
set(gca,'yscale','log');

%% 

[cfs,f] = cwt(x_2,eco_fs);
cwt(x_2,eco_fs);

%%
F1 = 19; F2 = 34;
cfs(f > F1 & f < F2) = 0;
xrec = icwt(cfs);
%Display the CWT of the reconstructed signal. The initial 25-Hz component is removed.
cwt(xrec)


