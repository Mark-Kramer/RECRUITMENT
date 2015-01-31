% Approximate the LFP as the summed synaptic activity over a square of size
% group size.
%
% INPUTS.
%  LL = the side of the data square.
%  synaptic = the synaptic dynamics.
%  group_size = the side of the group over which to average.
%
% OUTPUTS.
%  LFP = the approixmate LFP, of size [time, group_size, group_size].

function LFP = approximate_LFP(LL,synaptic,group_size)

  CI = zeros(LL,LL);

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
  LFP = zeros(size(synaptic,1),length(i0));
  for k=1:length(i0)
      ind0 = find(CI == i0(k));
      LFP(:,k) = mean(synaptic(:,ind0),2);
  end
  
end