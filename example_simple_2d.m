%%  Simulate RS cell in a 2-dim grid.

clear
close all

%L=2704; LL=52;
L=144; LL=12;           %L is the total # of neurons.
                        %LL is the size of a side of the 2-dim grid.

T = 5000;               %Total # of steps to take.

rs0.C = 1.0;            %HH model cell parameters.
rs0.gL  = 0.5;
rs0.gNaF = 150;
rs0.gKDR = 300;
rs0.gNaP = 0;
rs0.gKM  = 0;
rs0.gKNa = 0;
rs0.sigma = 1;
rs0.I0 = zeros(1,L);
rs0.I0(55) = -5;

syn0.taudEE=1.0;       %Synapse parameters.
syn0.taurEE=0.5;
syn0.gee=0.1;  

%Synaptic connections
p=0.0; [CEE, RCposition] = make_nn_network_full_square_small_world(L,LL,p);
CEE=CEE*1/8;                %Scale the connections so that total input from all nbrs is 1.
for k=1:L
    CEE(k,k)=0;             %No autapses.
end
C0.EtoE=CEE;

ic = 0;                 %No initial conditions to start.

%% Initial runs with synapses OFF.

syn0.gee=0;             %Synapses OFF.

for nsim=1:2            %Run a couple of times to initialize ...
    fprintf(['Initialize run ' num2str(nsim) '... \n'])
    [V,t,ic,current,synaptic] = simple_rs_2d(T,L,ic,  rs0,syn0,C0,RCposition);
end
fprintf('Initialization finished. \n')
imagesc(V, [-80 20]);
save('ic.mat', 'ic');   %Save the last step as the initial conditions, ar

%% Turn on connectivity.
syn0.gee=3;
load('ic.mat')
[V,t,ic,current,synaptic] = simple_rs_2d(T,L,ic,  rs0,syn0,C0,RCposition);
imagesc(V, [-80 20]);

%% Make connectivity small world.
p=0.01; [CEE, RCposition] = make_nn_network_full_square_small_world(L,LL,p);
CEE=CEE*1/8;                %Make long dist connections weaker.
for k=1:L
    CEE(k,k)=0;             %Make on module connections stronger.
end
C0.EtoE=CEE;
imagesc(C0.EtoE)

%% ... And run it.
syn0.gee=3;
load('ic.mat')
[V,t,ic,current,synaptic] = simple_rs_2d(T,L,ic,  rs0,syn0,C0,RCposition);
imagesc(V, [-80 20])

%% Make a movie.

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

%%  Extract average synaptic for a "group" of neurons.

CI = zeros(LL,LL);
group_size = 5;
counter=1;
for r=2:group_size:LL-1
    for c=2:group_size:LL-1
        CI(r:r+group_size-1,c:c+group_size-1) = counter;
        counter=counter+1;
    end
end

i0 = unique(CI);
good = find(i0 > 0);
i0 = i0(good);
LFP = zeros(size(V,1),length(i0));
for k=1:length(i0)
    ind0 = find(CI == i0(k));
    LFP(:,k) = mean(synaptic(:,ind0),2);
end

plot(LFP)