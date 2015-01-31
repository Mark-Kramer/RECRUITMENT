function [C, positionRC] = make_nn_network_full_square_small_world(LL,L,p)

  C = zeros(LL);

  %Top Left.
  i=1;
  j=1;                   %x,y  
  pre    = sub2ind([L,L], i,j);
  postL0 = sub2ind([L,L], L  ,j);
  postLU = sub2ind([L,L], L  ,j+1);
  postLD = sub2ind([L,L], L  ,L);
  postR0 = sub2ind([L,L], i+1,j);
  postRU = sub2ind([L,L], i+1,j+1);
  postRD = sub2ind([L,L], i+1,L);
  postU0 = sub2ind([L,L], i,j+1);
  postD0 = sub2ind([L,L], i,L);
  C(postL0, pre) = 1;
  C(postLU, pre) = 1;
  C(postLD, pre) = 1;
  C(postR0, pre) = 1;
  C(postRU, pre) = 1;
  C(postRD, pre) = 1;
  C(postU0, pre) = 1;
  C(postD0, pre) = 1;

  %Bottom Left.
  i=L;
  j=1;
  pre   = sub2ind([L,L], i,j);
  postL0 = sub2ind([L,L], i-1,j);
  postLU = sub2ind([L,L], i-1,j+1);
  postLD = sub2ind([L,L], i-1,L);  
  postR0 = sub2ind([L,L], 1,j);
  postRU = sub2ind([L,L], 1,j+1);
  postRD = sub2ind([L,L], 1,L);
  postU0 = sub2ind([L,L], i,j+1);
  postD0 = sub2ind([L,L], i,L);
  C(postL0, pre) = 1;
  C(postLU, pre) = 1;
  C(postLD, pre) = 1;
  C(postR0, pre) = 1;
  C(postRU, pre) = 1;
  C(postRD, pre) = 1;
  C(postU0, pre) = 1;
  C(postD0, pre) = 1;
  
  %Top Right.
  i=1;
  j=L;
  pre   = sub2ind([L,L], i,j);
  postL0 = sub2ind([L,L], L,j);
  postLU = sub2ind([L,L], L,1);
  postLD = sub2ind([L,L], L,j-1);
  postR0 = sub2ind([L,L], i+1,j);
  postRU = sub2ind([L,L], i+1,1);
  postRD = sub2ind([L,L], i+1,j-1);
  postU0 = sub2ind([L,L], i,1);
  postD0 = sub2ind([L,L], i,j-1);
  C(postL0, pre) = 1;
  C(postLU, pre) = 1;
  C(postLD, pre) = 1;
  C(postR0, pre) = 1;
  C(postRU, pre) = 1;
  C(postRD, pre) = 1;
  C(postU0, pre) = 1;
  C(postD0, pre) = 1;
  
  %Bottom Right.
  i=L;
  j=L;
  pre   = sub2ind([L,L], i,j);
  postL0 = sub2ind([L,L], i-1,j);
  postLU = sub2ind([L,L], i-1,1);
  postLD = sub2ind([L,L], i-1,j-1);
  postR0 = sub2ind([L,L], 1,j);
  postRU = sub2ind([L,L], 1,1);
  postRD = sub2ind([L,L], 1,j-1);
  postU0 = sub2ind([L,L], i,1);
  postD0 = sub2ind([L,L], i,j-1);
  C(postL0, pre) = 1;
  C(postLU, pre) = 1;
  C(postLD, pre) = 1;
  C(postR0, pre) = 1;
  C(postRU, pre) = 1;
  C(postRD, pre) = 1;
  C(postU0, pre) = 1;
  C(postD0, pre) = 1;

  %Top Row.
  i=1;
  for j=2:L-1
      pre   = sub2ind([L,L], i,j);
      postL0 = sub2ind([L,L], L,j);
      postLU = sub2ind([L,L], L,j+1);
      postLD = sub2ind([L,L], L,j-1);
      postR0 = sub2ind([L,L], i+1,j);
      postRU = sub2ind([L,L], i+1,j+1);
      postRD = sub2ind([L,L], i+1,j-1);
      postU0 = sub2ind([L,L], i,j+1);
      postD0 = sub2ind([L,L], i,j-1);
      C(postL0, pre) = 1;
      C(postLU, pre) = 1;
      C(postLD, pre) = 1;
      C(postR0, pre) = 1;
      C(postRU, pre) = 1;
      C(postRD, pre) = 1;
      C(postU0, pre) = 1;
      C(postD0, pre) = 1;
  end
  
  %Bottom Row.
  i=L;
  for j=2:L-1
      pre   = sub2ind([L,L], i,j);
      postL0 = sub2ind([L,L], i-1,j);
      postLU = sub2ind([L,L], i-1,j+1);
      postLD = sub2ind([L,L], i-1,j-1);
      postR0 = sub2ind([L,L], 1,j);
      postRU = sub2ind([L,L], 1,j+1);
      postRD = sub2ind([L,L], 1,j-1);
      postU0 = sub2ind([L,L], i,j+1);
      postD0 = sub2ind([L,L], i,j-1);
      C(postL0, pre) = 1;
      C(postLU, pre) = 1;
      C(postLD, pre) = 1;
      C(postR0, pre) = 1;
      C(postRU, pre) = 1;
      C(postRD, pre) = 1;
      C(postU0, pre) = 1;
      C(postD0, pre) = 1;
  end
  
  %Left Row.
  j=1;
  for i=2:L-1
      pre   = sub2ind([L,L], i,j);
      postL0 = sub2ind([L,L], i-1,j);
      postLU = sub2ind([L,L], i-1,j+1);
      postLD = sub2ind([L,L], i-1,L);
      postR0 = sub2ind([L,L], i+1,j);
      postRU = sub2ind([L,L], i+1,j+1);
      postRD = sub2ind([L,L], i+1,L);
      postU0 = sub2ind([L,L], i,j+1);
      postD0 = sub2ind([L,L], i,L);
      C(postL0, pre) = 1;
      C(postLU, pre) = 1;
      C(postLD, pre) = 1;
      C(postR0, pre) = 1;
      C(postRU, pre) = 1;
      C(postRD, pre) = 1;
      C(postU0, pre) = 1;
      C(postD0, pre) = 1;
  end
  
  %Right Row.
  j=L;
  for i=2:L-1
      pre   = sub2ind([L,L], i,j);
      postL0 = sub2ind([L,L], i-1,j);
      postLU = sub2ind([L,L], i-1,1);
      postLD = sub2ind([L,L], i-1,j-1);
      postR0 = sub2ind([L,L], i+1,j);
      postRU = sub2ind([L,L], i+1,1);
      postRD = sub2ind([L,L], i+1,j-1);
      postU0 = sub2ind([L,L], i,1);
      postD0 = sub2ind([L,L], i,j-1);
      C(postL0, pre) = 1;
      C(postLU, pre) = 1;
      C(postLD, pre) = 1;
      C(postR0, pre) = 1;
      C(postRU, pre) = 1;
      C(postRD, pre) = 1;
      C(postU0, pre) = 1;
      C(postD0, pre) = 1;
  end
  
  %Middle parts.
  for i=2:L-1
      for j=2:L-1
          pre   = sub2ind([L,L], i,j);
          postL0 = sub2ind([L,L], i-1,j);
          postLU = sub2ind([L,L], i-1,j+1);
          postLD = sub2ind([L,L], i-1,j-1);
          postR0 = sub2ind([L,L], i+1,j);
          postRU = sub2ind([L,L], i+1,j+1);
          postRD = sub2ind([L,L], i+1,j-1);
          postU0 = sub2ind([L,L], i,j+1);
          postD0 = sub2ind([L,L], i,j-1);
          
          C(postL0, pre) = 1;
          C(postLU, pre) = 1;
          C(postLD, pre) = 1;
          C(postR0, pre) = 1;
          C(postRU, pre) = 1;
          C(postRD, pre) = 1;
          C(postU0, pre) = 1;
          C(postD0, pre) = 1;
      end
  end

  counter=1;
  for i=1:L
      for j=1:L
          positionRC(counter,1)=i;
          positionRC(counter,2)=j;
          counter=counter+1;
      end
  end
  
  %Make small world.
  %Find each (i,j) pair where i>j --- that's all node pairs!
  [row, col] = find(triu(C));
  for k=1:length(row)           %For each (i,j) pair.
      i=row(k);
      j=col(k);
      if rand() < p             %Draw a random number [0,1].  If it's < p,
          C(i,j)=0;             %  then eliminate this edge (i,j),
          C(j,i)=0;             %  and eliminate (j,i).
          inew=1;
          jnew=1;               %Now find a new edge.
          while inew == jnew || C(inew,jnew) ==1        %If i==j or the edge already exists,
              ijnew = randperm(LL);                      %Then guess a new (i,j) pair.
              inew = ijnew(1);                          %Repeat until we find a new possible edge.
              jnew = ijnew(2);
          end
          C(inew,jnew)=1;       %Put the new edge in the network.
          C(jnew,inew)=1;
      end
  end
  
end
