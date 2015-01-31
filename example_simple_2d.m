%%  Simulate RS and FS.

clear
close all

%L=2704; LL=52;
L=100; LL=10;

T = 5000;

rs0.C = 1.0;
rs0.gL  = 0.5;
rs0.gNaF = 150;
rs0.gKDR = 300;
rs0.gNaP = 0.0;
rs0.gKM  = 0;
rs0.gKNa = 0;
rs0.sigma = 1;
rs0.I0 = zeros(1,L);
rs0.I0(55) = -5;

syn0.taudEE=1.0;       %10.0;
syn0.taurEE=0.5;

syn0.tauGEE=0;              %Turn OFF facilitation.

syn0.gee=0.1;                 %Facilitates, so be careful.

%Structural connections
%m[CEE, RCposition] = make_nn_network_full_square(L,LL);
p=0.0; [CEE, RCposition] = make_nn_network_full_square_small_world(L,LL,p);

CEE=CEE*1/8;                %Make long dist connections weaker.
for k=1:L
    CEE(k,k)=0;             %Make on module connections stronger.
end

% %Within-module connections only.
% Cin   = zeros(L);
% for k=1:L
%     Cin(k,k)=1;
% end

C0.EtoE=CEE;
%C0.ItoE=Cin;
%C0.ItoI=Cin;

ic=0;
V_rs=[];


%% RS only, initial runs.
EK=-80;
syn0.gee=0;
[V,t,ic,current,synaptic] = simple_rs_2d(T,L,EK,ic,  rs0,syn0,C0,RCposition);
imagesc(V, [-80 20]);
save('ic.mat', 'ic');

%% RS only, turn on connectivity.
syn0.gee=3;
load('ic.mat')
[V,t,ic,current,synaptic] = simple_rs_2d(T,L,EK,ic,  rs0,syn0,C0,RCposition);
imagesc(V, [-80 20]);

%% Make small world.
p=0.01; [CEE, RCposition] = make_nn_network_full_square_small_world(L,LL,p);
CEE=CEE*1/8;                %Make long dist connections weaker.
for k=1:L
    CEE(k,k)=0;             %Make on module connections stronger.
end
C0.EtoE=CEE;
imagesc(C0.EtoE)

%% And run it.
syn0.gee=3;
load('ic.mat')
[V,t,ic,current,synaptic] = simple_rs_2d(T,L,EK,ic,  rs0,syn0,C0,RCposition);
imagesc(V, [-80 20])

%%  Movie.

i0=1;
Lside = sqrt(size(V,2));
for i=i0:100:size(V,1)
    map = reshape(squeeze(V(i,:)), [Lside,Lside]);
    map = map(2:Lside-1, 2:Lside-1);
    imagesc(map, [-80, 40])
    title(num2str(t(i)))
    pause(0.1)
    %keyboard
end

%%  Extract average synaptic.

i0=1;
Lside = sqrt(size(V,2));
%for i=i0:100:size(V,1)
i=100;
     map = reshape(squeeze(synaptic(i,:)), [Lside,Lside]);
     map = map(2:Lside-1, 2:Lside-1);

%%
CI = zeros(LL,LL);
group_size = 5;
counter=1;
for r=1:group_size:LL
    for c=1:group_size:LL
        CI(r:r+group_size-1,c:c+group_size-1) = counter;
        counter=counter+1;
    end
end



%%  FS and RS cell (testing).

% %syn0.gei=10.0;
% %ic.GEE=5.0;
% EK0=ones(1,L)*-80;
% [V0,rs_V,t,mNaF,hNaF,mKDR,mCaH,kV,mKM,ic] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
% % V = [V;V0];
% % V_rs = [V_rs; rs_V];
% % I = [I,I0];
% % EK= [EK,EK0];
% 
% subplot(2,1,1)
% plot(V0)
% ylim([-100 30])
% subplot(2,1,2)
% plot(rs_V)
% ylim([-100,30])

%%  FS and RS run it.

% %Load the ICs.
% load('/Users/mak/research/OMAR/dat/good_ic.mat')
% 
% 
% syn0.taudEE=10.0;       %10.0
% syn0.taurEE=5.0;        %5.0
% syn0.taudEI=10.0;
% syn0.taurEI=0.5;
% 
% syn0.taudII=10.0;
% syn0.taudIE=10.0;
% 
% syn0.tauGEE=0;              %Turn OFF facilitation.
% 
% syn0.gee=10;                 %Facilitates, so be careful.
% syn0.gei=0.25;
% syn0.gie=1.25;
% syn0.gii=0;

EK0=ones(1,L)*-55; 
ic.GEE=ones(1,L)*100;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);

imagesc(rs_V)

%%


EK0=ones(1,L)*-90;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
ic.GEE=ones(1,L)*1;
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);
save(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' num2str(EK0(1)) '_' num2str(ic.GEE(1)) '.mat'], 'fs_V', 'rs_V', 'LFP', 't', 'EK0', 'ic')
clear('fs_V', 'rs_V', 'LFP', 't');

EK0=ones(1,L)*-80;
ic.GEE=ones(1,L)*1;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);
save(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' num2str(EK0(1)) '_' num2str(ic.GEE(1)) '.mat'], 'fs_V', 'rs_V', 'LFP', 't', 'EK0', 'ic')
clear('fs_V', 'rs_V', 'LFP', 't');

EK0=ones(1,L)*-70;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
ic.GEE=ones(1,L)*1;
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);
save(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' num2str(EK0(1)) '_' num2str(ic.GEE(1)) '.mat'], 'fs_V', 'rs_V', 'LFP', 't', 'EK0', 'ic')
clear('fs_V', 'rs_V', 'LFP', 't');

EK0=ones(1,L)*-60;
ic.GEE=ones(1,L)*1;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);
save(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' num2str(EK0(1)) '_' num2str(ic.GEE(1)) '.mat'], 'fs_V', 'rs_V', 'LFP', 't', 'EK0', 'ic')
clear('fs_V', 'rs_V', 'LFP', 't');

EK0=ones(1,L)*-55;
ic.GEE=ones(1,L)*1;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);
save(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' num2str(EK0(1)) '_' num2str(ic.GEE(1)) '.mat'], 'fs_V', 'rs_V', 'LFP', 't', 'EK0', 'ic')
clear('fs_V', 'rs_V', 'LFP', 't');

EK0=ones(1,L)*-55;
ic.GEE=ones(1,L)*10;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);
save(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' num2str(EK0(1)) '_' num2str(ic.GEE(1)) '.mat'], 'fs_V', 'rs_V', 'LFP', 't', 'EK0', 'ic')
clear('fs_V', 'rs_V', 'LFP', 't');

EK0=ones(1,L)*-55;
ic.GEE=ones(1,L)*25;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);
save(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' num2str(EK0(1)) '_' num2str(ic.GEE(1)) '.mat'], 'fs_V', 'rs_V', 'LFP', 't', 'EK0', 'ic')
clear('fs_V', 'rs_V', 'LFP', 't');

EK0=ones(1,L)*-55;
ic.GEE=ones(1,L)*50;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);
save(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' num2str(EK0(1)) '_' num2str(ic.GEE(1)) '.mat'], 'fs_V', 'rs_V', 'LFP', 't', 'EK0', 'ic')
clear('fs_V', 'rs_V', 'LFP', 't');

EK0=ones(1,L)*-55;
ic.GEE=ones(1,L)*75;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);
save(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' num2str(EK0(1)) '_' num2str(ic.GEE(1)) '.mat'], 'fs_V', 'rs_V', 'LFP', 't', 'EK0', 'ic')
clear('fs_V', 'rs_V', 'LFP', 't');

EK0=ones(1,L)*-55;
ic.GEE=ones(1,L)*100;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);
save(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' num2str(EK0(1)) '_' num2str(ic.GEE(1)) '.mat'], 'fs_V', 'rs_V', 'LFP', 't', 'EK0', 'ic')
clear('fs_V', 'rs_V', 'LFP', 't');

EK0=ones(1,L)*-55;
ic.GEE=ones(1,L)*125;
[fs_V,rs_V,t,~,~,~,~,~,~,ic,~,LFP] = traub_fs_rs_v5_2d(T,L, I0, gLNa, gLK, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK0, C, sigma,ic, rs0,syn0,C0,RCposition);
fs_V = downsample(fs_V,10);
rs_V = downsample(rs_V,10);
LFP  = downsample(LFP,10);
t = downsample(t,10);
save(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' num2str(EK0(1)) '_' num2str(ic.GEE(1)) '.mat'], 'fs_V', 'rs_V', 'LFP', 't', 'EK0', 'ic')
clear('fs_V', 'rs_V', 'LFP', 't');

%%  Movie.

i0=1;
Lside=10;
V = LFP;
Lside = sqrt(size(V,2));
for i=i0:10:size(V,1)
    map = reshape(squeeze(V(i,:)), [Lside,Lside]);
    map = map(2:Lside-1, 2:Lside-1);
    imagesc(map, [-400, 900])
    title(num2str(t(i)))
    pause(0.1)
    %keyboard
end

%%  Combined activity
EKnames = {'-100'; '-90'; '-80'; '-70'; '-60'; '-55'; '-55';'-55';'-55';'-55'; '-55'; '-55'};
GEnames = {   '1';   '1';   '1';   '1';   '1';   '1';  '10'; '25'; '50'; '75'; '100'; '125'};

rs_spks = zeros(length(EKnames)*10000,100);
fs_spks = zeros(length(EKnames)*10000,100);

LFP0 = [];
for nms=1:length(EKnames)
    fprintf([num2str(nms) '\n'])
    load(['~/research/OMAR/dat/traub_fs_rs_v5_2d_' EKnames{nms} '_' GEnames{nms} '.mat'])
    LFP0 = [LFP0; mean(LFP,2)];
    
    for k=1:size(rs_V,2);
        [~,i0] = findpeaks(rs_V(:,k), 'minpeakheight', 0);
        rs_spks(i0+(nms-1)*10000,k)=1;
    end
    
    for k=1:size(fs_V,2);
        [~,i0] = findpeaks(fs_V(:,k), 'minpeakheight', 0);
        fs_spks(i0+(nms-1)*10000,k)=1;
    end
    
end

%%
dt = t(10)-t(9);
taxis = (1:length(LFP0))*dt;

subplot(3,1,1)
plot(taxis,mean(LFP0,2))
subplot(3,1,2)
rs_fr = mean(rs_spks,2)/dt;
plot(taxis, smooth(rs_fr,101))
fs_fr = mean(fs_spks,2)/dt;
hold on
plot(taxis, smooth(fs_fr,101), 'g')
hold off

dT=500;
tbinned = (1:dT:max(taxis)-dT);

for i=1:length(tbinned)
    igood = find(taxis > tbinned(i) & taxis < tbinned(i)+dT);
    rs_n(i) = sum(sum(rs_spks(igood,:)))/size(rs_V,2);
    fs_n(i) = sum(sum(fs_spks(igood,:)))/size(fs_V,2);
end
subplot(3,1,3)
plot(tbinned+dT/2, rs_n, 'LineWidth', 4)
hold on
plot(tbinned+dT/2, fs_n, 'g', 'LineWidth', 4)
hold off


