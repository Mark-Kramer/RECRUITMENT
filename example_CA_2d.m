%% Simulate cellular automata in a 2-dim grid.

clear
%close all

%L is the total # of neurons.
%LL is the size of a side of the 2-dim grid.
%For example, with LL=27, and L=27^2=729, we'll eventually
%create a 25-by-25 grid, and drop the outer boundary.
L=144; LL=12;

%Main parameters.
%p = probability of a rewiring a neighborly connection and making it a long
%distance connection.  When p=0, the network is neighbor-to-neighbor.  When
%p=1, the network is randomly connected.
%
%gee = overall strength of synapse.  The larger this value, the more likely
%that the input to a cell will "activate" it (i.e., flip it from 0 to 1).
p=0.1;
g=0.05;

%Make synaptic connections
[C, RCposition] = make_nn_network_full_square_small_world(L,LL,p);
C=C*1/8;                %Scale the connections so that total input from all nbrs is 1.
for k=1:L
    C(k,k)=0;           %No autapses.
end
C=g*C;
C0.EtoE=C;

initial_on_index = 14;      %Turn "on" an initial cell at this index.
t = simple_CA_2d(L,initial_on_index,C0,RCposition);
fprintf(['Took ' num2str(max(t)) ' steps to recruit entire grid \n'])

% plot(t,'r')
% hold on
% figure
%% Make a movie.

CA = zeros(LL,LL);
CA(initial_on_index) = 1;
map = CA(2:LL-1, 2:LL-1);
imagesc(map, [0,1])
pause(0.5)
for i=1:max(t)
    i0 = find(t==i);
    CA(i0)=1;
    map = CA(2:LL-1, 2:LL-1);
    imagesc(map, [0,1])
    title(num2str(i))
    pause(0.001)
    %keyboard
end


