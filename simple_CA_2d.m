function [t] = simple_CA_2d(L,initial_on_index,C0,positionRC)
  
  CA = zeros(1,L);
  CA(initial_on_index) = 1;
  
  % Prepare for no flux.
  row1 = find(positionRC(:,1)==1);
  row2 = find(positionRC(:,1)==2);
  lastRow = max(positionRC(:,1));
  rowEnd   = find(positionRC(:,1)==lastRow);
  rowEndm1 = find(positionRC(:,1)==lastRow-1);
  
  col1 = find(positionRC(:,2)==1);
  col2 = find(positionRC(:,2)==2);
  lastCol = max(positionRC(:,2));
  colEnd   = find(positionRC(:,2)==lastCol);
  colEndm1 = find(positionRC(:,2)==lastCol-1);
      
  step=0;
  t = nan(1,L);
  t(initial_on_index)=step;
  step=step+1;
  while sum(CA) < L
      i0 = find(CA == 0);
      input = CA * C0.EtoE;
      CA0 = binornd(1,input,1,L);
      CA(i0) = CA0(i0);
      ichanged = i0(find(CA0(i0)));
      t(ichanged)=step;
      step=step+1;
         
      %No flux.      
      CA(row1) = CA(row2);      
      CA(rowEnd) = CA(rowEndm1);
      CA(col1) = CA(col2);
      CA(colEnd) = CA(colEndm1);
      
      %plot(CA)
  end
      

      
end    
