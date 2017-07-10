function main_CA_2d

LL=10; L=LL^2;
nbSz = 4;
nbSeed = 32;
% nbSz = 1;
% nbSeed = 1;
randSeed = false;
g = 0.05;
pr = 0;
psw = 0.08;
figName = 'sim';
global MAKEMOVIE;
MAKEMOVIE = false;
global SIMFOLDER;
SIMFOLDER = 'toolboxes/mark/model/sim/sim_step';

save('paramvars.mat');
params = load('paramvars.mat'); %#ok<NASGU>
delete('paramvars.mat');

%% testing the connectivity param to be in the SW range
% p = [0 0.0001:0.0002:0.001 0.002:0.002:0.01 0.02:0.02:0.1 0.2:0.2:1];
% LEN = zeros(1,length(p));
% COPEN = zeros(1,length(p));
% len = zeros(1,20);
% copen = zeros(1,20);
% for i = 1:length(p)
%     for j = 1 : 20
%         C = make_nn_network_full_square_small_world(L,LL,p(i));
%         [len(j),~,~,~,copen(j),~] = graphProperties(C);
%     end
%     len = len(isfinite(len));
%     copen = copen(isfinite(copen));
%     LEN(i) = mean(len);
%     COPEN(i) = mean(copen);
%     fprintf('%d\n',i);
% end
% figure; plotyy(p,LEN/LEN(1),p,COPEN/COPEN(1),@semilogx);
% saveas(gcf, 'sw_param', 'pdf');


%%
idx = reshape(1:L,LL,LL);
idx = idx(2:end-1,2:end-1);
if randSeed
    % random seed
    initIdx = idx(:);
    initIdx = idx(randi(length(idx), 1, nbSeed));
else
    % distributed seed
    initIdx = idx(checkerboard(1,LL/2-1,LL/2-1)>0.5);
end
% figure; M = zeros(LL); M(initIdx) = 1; imagesc(M);

Reg.T = cell(1,nbSeed);
Reg.MORAN = cell(1,nbSeed);
Reg.CORR = cell(1,nbSeed);
for i = 1 : nbSeed
    [Reg.T{i}, Reg.MORAN{i}, Reg.CORR{i}] = runsimulation(LL,L,pr,g,initIdx(i),nbSz);
end

SW.T = cell(1,nbSeed);
SW.MORAN = cell(1,nbSeed);
SW.CORR = cell(1,nbSeed);
for i = 1 : nbSeed
    [SW.T{i}, SW.MORAN{i}, SW.CORR{i}] = runsimulation(LL,L,psw,g,initIdx(i),nbSz);
end
save([figName '.mat'], 'Reg', 'SW', 'params');
% return;

%%
figName = [figName '.pdf'];

m = cellfun(@mean, SW.MORAN);
mr = cellfun(@mean, Reg.MORAN);
plotbarstats(mr, m, 'Moran''s Index', figName);

c = cellfun(@mean, SW.CORR);
cr = cellfun(@mean, Reg.CORR);
plotbarstats(cr, c, 'Correlation between maps', figName);

t = cellfun(@mean, SW.T);
tr = cellfun(@mean, Reg.T);
plotbarstats(tr, t, 'Recruitment time', figName);


end

function plotbarstats(grp1, grp2, ytitle, figName)

f = figure;

bar([mean(grp1) mean(grp2)]);
hold on; 
errorbar([mean(grp1) mean(grp2)], [std(grp1) std(grp2)], '.');
set(gca, 'xtickLabel', {'Regular', 'Small-world'})
xlim([0 3])
ylabel(ytitle)
[p,tab]=anova1([grp1,grp2], [ones(1,length(grp1)) 2*ones(1,length(grp2))], 'off');
title(['F(' num2str(tab{2,3}) ',' num2str(tab{3,3}) ')=' num2str(tab{2,5}) ' / p=' num2str(p)]);

export_fig('-append', figName, f);
close(f);

end


function [T, MORAN, CORR] = runsimulation(LL,L,p,g,initial_on_index, nbSz)

% Simulate cellular automata in a 2-dim grid.
%
% L is the total # of neurons.
% LL is the size of a side of the 2-dim grid.
%For example, with LL=27, and L=27^2=729, we'll eventually
%create a 25-by-25 grid, and drop the outer boundary.
%
%p = probability of a rewiring a neighborly connection and making it a long
%distance connection.  When p=0, the network is neighbor-to-neighbor.  When
%p=1, the network is randomly connected.
%
%g = overall strength of synapse.  The larger this value, the more likely
%that the input to a cell will "activate" it (i.e., flip it from 0 to 1).
% 
% initial_on_index = 14;      %Turn "on" an initial cell at this index.

%Make synaptic connections
[C, RCposition] = make_nn_network_full_square_small_world(L,LL,p);

global MAKEMOVIE;
global SIMFOLDER;

C=C*1/8;                %Scale the connections so that total input from all nbrs is 1.
for k=1:L
    C(k,k)=0;           %No autapses.
end
C=g*C;
C0.EtoE=C;

MAP = zeros(L-2*LL-2*(LL-2),nbSz);
T = zeros(1,nbSz);
MORAN = zeros(1,nbSz);
for i = 1:nbSz
    t = simple_CA_2d(L,initial_on_index,C0,RCposition);
    T(i) = max(t);
    map = reshape(t,LL,LL);
    map = map(2:LL-1, 2:LL-1);
    moran = localmoran(map, 'neighborCircle', 3);
	MORAN(i) = nansum(moran(:));
    MAP(:,i) = map(:);
    % % Plots
    % fprintf(['Took ' num2str(max(t)) ' steps to recruit entire grid \n'])
    % plot(t,'r')
    % figure, imagesc(-map);

    % Make a movie.
    if MAKEMOVIE
        f = figure;
        CA = zeros(LL,LL);
        CA(initial_on_index) = 1;
        map = CA(2:LL-1, 2:LL-1);
        imagesc(map, [0,1])
    %     pause(0.5)
        saveas(f, [SIMFOLDER 'sim_step0'], 'png');
        for step=1:max(t)
            CA(t==step)=1;
            map = CA(2:LL-1, 2:LL-1);
            imagesc(map, [0,1])
            colormap(1-gray);
    %         title(num2str(i))
    %         pause(0.001)
            %keyboard
            axis off
            saveas(f, [SIMFOLDER num2str(step)], 'png');
        end
        close(f)
    end
end
CORR = corr(MAP, 'type', 'pearson', 'rows', 'pairwise');
CORR = CORR(triu(true(size(CORR)),1));

end
