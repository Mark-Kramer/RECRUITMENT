%NOTES.
%  From [Cunningham PNAS 2004 SI].
%  x.  No persistent Na
%  x.  No h-current
%  x.  Updated fast Na activation (done)
%  x.  Updated fast Na inactivation (done)
%  x.  Updated KDR (done)
%
%April 3, 2014.  Fixed Cunningham 2004 typo in alpha_mNaF, beta_mNaF
%  function.
%
%April 3, 2014.  There's a problem with the "noise" term in FS cell.  When
%  the noise is just above 0, the FS cell gets fat spikes.  

function [rs_V,t,ic,current,synaptic] = simple_WC_2d(T,L,ic,  rs0,syn0,C0,positionRC)

  EK = -80;

  dt = 0.005;
  
  rs_V = zeros(T,L);
  rs_mNaF = zeros(T,L);
  rs_hNaF = zeros(T,L);
  rs_mNaP = zeros(T,L);
  rs_mKDR = zeros(T,L);
  rs_mKM = zeros(T,L);
  rs_nKNa = zeros(T,L);
   
  sEE = zeros(T,L);
  %GEE = zeros(T,L);  GEE(1,:)=syn0.gee;
  GEE = syn0.gee;
  
  current = zeros(T,L);
  synaptic= zeros(T,L);
  
  if isstruct(ic)      
      rs_V(1,:) = ic.rs_V;
      rs_mNaF(1,:) = ic.rs_mNaF;
      rs_mNaP(1,:) = ic.rs_mNaP;
      rs_hNaF(1,:) = ic.rs_hNaF;
      rs_mKDR(1,:) = ic.rs_KDR;
      rs_mKM(1,:)  = ic.rs_mKM;
      rs_nKNa(1,:) = ic.rs_nKNa;
             
      sEE(1,:) = ic.sEE;
      
      %GEE(1,:) = ic.GEE;
      current(1,:) = ic.current;
  end
  
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
      
  for i=1:T-1
      
      %       %.*.*.*.*  Synapses.  .*.*.*.*;
      sEE(i+1,:) = sEE(i,:) + dt.*(-sEE(i,:)./syn0.taudEE + (1 - sEE(i,:))   .*(1 + tanh(rs_V(i,:)./10.))./syn0.taurEE);
      sEtoE = sEE(i,:) * C0.EtoE;
      %GEE(i+1,:) = min(GEE(i,:) + dt*syn0.tauGEE*(0.01*(0-GEE(i,:))+3*sEtoE), 100);
      %GEE(i+1,:) = GEE(i,:) + dt*(syn0.tauGEEd*(0-GEE(i,:))+syn0.tauGEEr*sEtoE);
  
      %%%% %%%% RS %%%% %%%%
      noise   = sqrt(dt).*randn(1,L).*rs0.sigma;

      rs_V(i+1,:) = rs_V(i,:) + dt.*rs0.C*( -rs0.I0  ...
          - rs0.gL*(70.0 + rs_V(i,:)) ...
          - (  rs0.gNaF.*rs_mNaF(i,:).^3.0.*rs_hNaF(i,:) + ...
               rs0.gNaP.*rs_mNaP(i,:)).*(-50.0 + + rs_V(i,:)) ...             ;NaF current.
          - (  rs0.gKDR.*rs_mKDR(i,:).^4.0 ...                            ;KDR current.
             + rs0.gKM .*rs_mKM(i,:) ...                                  ;KDR current.
                                  ).*(-EK + rs_V(i,:)) ...   ;KNa current (Slack-B, slow)
             - GEE .*sEtoE.*    rs_V(i,:)) ...
             ... - syn0.gee *sEtoE.*(   rs_V(i,:)) ...
             ... - syn0.gie *sItoE.*(80+rs_V(i,:))) ...
             + noise;
         
      current(i,:) = (-rs0.I0  ...
          - rs0.gL*(70.0 + rs_V(i,:)) ...
          - (  rs0.gNaF.*rs_mNaF(i,:).^3.0.*rs_hNaF(i,:) + ...
               rs0.gNaP.*rs_mNaP(i,:)).*(-50.0 + + rs_V(i,:)) ...             ;NaF current.
          - (  rs0.gKDR.*rs_mKDR(i,:).^4.0 ...                            ;KDR current.
             + rs0.gKM .*rs_mKM(i,:) ...                                  ;KDR current.
                                  ).*(-EK + rs_V(i,:)) ...   ;KNa current (Slack-B, slow)
             - GEE  .*sEtoE.*    rs_V(i,:));  % ...
             %- syn0.gie *sItoE.*(80+rs_V(i,:)));
         
      synaptic(i,:) = ( ...
             - GEE .*sEtoE.*    rs_V(i,:));  % ...
             %- syn0.gie *sItoE.*(80+rs_V(i,:)));
      
      rs_mNaF(i+1,:) = rs_mNaF(i,:) + dt*(rs_alpha_m(rs_V(i,:)).*(1-rs_mNaF(i,:))  - rs_beta_m(rs_V(i,:)).*rs_mNaF(i,:));
      rs_hNaF(i+1,:) = rs_hNaF(i,:) + dt*(rs_alpha_h(rs_V(i,:)).*(1-rs_hNaF(i,:))  - rs_beta_h(rs_V(i,:)).*rs_hNaF(i,:));
      rs_mNaP(i+1,:) = rs_mNaP(i,:) + dt*(rs_alpha_mNaP(rs_V(i,:)).*(1-rs_mNaP(i,:))  - rs_beta_mNaP(rs_V(i,:)).*rs_mNaP(i,:));
      rs_mKDR(i+1,:) = rs_mKDR(i,:) + dt*(rs_alpha_n(rs_V(i,:)).*(1-rs_mKDR(i,:))  - rs_beta_n(rs_V(i,:)).*rs_mKDR(i,:));
      rs_mKM(i+1,:)  = rs_mKM(i,:)  + dt*(alpha_mKM( rs_V(i,:)).*(1- rs_mKM(i,:))  - beta_mKM( rs_V(i,:)).*rs_mKM(i,:));
      rs_nKNa(i+1,:) = rs_nKNa(i,:) + dt*(rs_alpha_nKNa(rs_V(i,:)).*(1-rs_nKNa(i,:)) - rs_beta_nKNa(rs_V(i,:)).*rs_nKNa(i,:));
     
      %No flux.      
      rs_V(i+1,row1) = rs_V(i+1,row2);
      rs_mNaF(i+1,row1) = rs_mNaF(i+1,row2);
      rs_hNaF(i+1,row1) = rs_hNaF(i+1,row2);
      rs_mNaP(i+1,row1) = rs_mNaP(i+1,row2);
      rs_mKDR(i+1,row1) = rs_mKDR(i+1,row2);
      rs_mKM(i+1,row1)  = rs_mKM(i+1,row2);
      rs_nKNa(i+1,row1) = rs_nKNa(i+1,row2);
      
      sEE(i+1,row1) = sEE(i+1,row2);
      current(i,row1) = current(i,row2);
      synaptic(i,row1)= synaptic(i,row2);

      rs_V(i+1,rowEnd) = rs_V(i+1,rowEndm1);
      rs_mNaF(i+1,rowEnd) = rs_mNaF(i+1,rowEndm1);
      rs_hNaF(i+1,rowEnd) = rs_hNaF(i+1,rowEndm1);
      rs_mNaP(i+1,rowEnd) = rs_mNaP(i+1,rowEndm1);
      rs_mKDR(i+1,rowEnd) = rs_mKDR(i+1,rowEndm1);
      rs_mKM(i+1,rowEnd)  = rs_mKM(i+1,rowEndm1);
      rs_nKNa(i+1,rowEnd) = rs_nKNa(i+1,rowEndm1);

      sEE(i+1,rowEnd) = sEE(i+1,rowEndm1);
      current(i,rowEnd) = current(i,rowEndm1);
      synaptic(i,rowEnd)= synaptic(i,rowEndm1);
      
      rs_V(i+1,col1) = rs_V(i+1,col2);
      rs_mNaF(i+1,col1) = rs_mNaF(i+1,col2);
      rs_hNaF(i+1,col1) = rs_hNaF(i+1,col2);
      rs_mNaP(i+1,col1) = rs_mNaP(i+1,col2);
      rs_mKDR(i+1,col1) = rs_mKDR(i+1,col2);
      rs_mKM(i+1,col1)  = rs_mKM(i+1,col2);
      rs_nKNa(i+1,col1) = rs_nKNa(i+1,col2);
   
      sEE(i+1,col1) = sEE(i+1,col2);
      %GEE(i+1,col1) = GEE(i+1,col2);
      current(i,col1) = current(i,col2);
      synaptic(i,col1)= synaptic(i,col2);
      
      rs_V(i+1,colEnd) = rs_V(i+1,colEndm1);
      rs_mNaF(i+1,colEnd) = rs_mNaF(i+1,colEndm1);
      rs_hNaF(i+1,colEnd) = rs_hNaF(i+1,colEndm1);
      rs_mNaP(i+1,colEnd) = rs_mNaP(i+1,colEndm1);
      rs_mKDR(i+1,colEnd) = rs_mKDR(i+1,colEndm1);
      rs_mKM(i+1,colEnd)  = rs_mKM(i+1,colEndm1);
      rs_nKNa(i+1,colEnd) = rs_nKNa(i+1,colEndm1);

      sEE(i+1,colEnd) = sEE(i+1,colEndm1);
      %GEE(i+1,colEnd) = GEE(i+1,colEndm1);
      current(i,colEnd) = current(i,colEndm1);
      synaptic(i,colEnd)= synaptic(i,colEndm1);
      
  end
  
  %Last current step.
  current(i+1,:) = (-rs0.I0  ...
          - rs0.gL*(70.0 + rs_V(i+1,:)) ...
          - (  rs0.gNaF.*rs_mNaF(i+1,:).^3.0.*rs_hNaF(i+1,:) + ...
               rs0.gNaP.*rs_mNaP(i+1,:)).*(-50.0 + + rs_V(i+1,:)) ...             ;NaF current.
          - (  rs0.gKDR.*rs_mKDR(i+1,:).^4.0 ...                            ;KDR current.
             + rs0.gKM .*rs_mKM(i+1,:) ...                                  ;KDR current.
                                  ).*(-EK + rs_V(i+1,:)) ...   ;KNa current (Slack-B, slow)
             - GEE .*sEtoE.*    rs_V(i+1,:)); %...
             %- syn0.gie *sItoE.*(80+rs_V(i+1,:)));
  current(i+1,row1) = current(i+1,row2);
  current(i+1,rowEnd) = current(i+1,rowEndm1);
  current(i+1,col1) = current(i+1,col2);
  current(i+1,colEnd) = current(i+1,colEndm1);
  
  
  ic = {};
  t = dt*(1:T);
  
  ic.rs_V = rs_V(end,:);
  ic.rs_mNaF = rs_mNaF(end,:);
  ic.rs_hNaF = rs_hNaF(end,:);
  ic.rs_mNaP = rs_mNaP(end,:);
  ic.rs_KDR  = rs_mKDR(end,:);
  ic.rs_mKM  = rs_mKM(end,:);
  ic.rs_nKNa = rs_nKNa(end,:);

  ic.sEE = sEE(end,:);
  %ic.GEE = GEE(end,:);
  ic.current = current(end,:);

end
    
function res = alpha_mNaF(V)    % NaF activation [Cunningham SI 2004].
  minf = 1.0 ./ (1.0 + exp( (-V-38.0)/10 ));
  taum = 0.0163 + 0.1487 * exp(-abs(-V-30.0)/10);
  res = minf./taum;
end

function res = beta_mNaF(V)    % NaF activation  [Cunningham SI 2004].
  minf = 1.0 ./ (1.0 + exp( (-V-38.0)/10 ));
  taum = 0.0163 + 0.1487 * exp(-abs(-V-30.0)/10);
  res = (1.0 - minf)./taum;
end

function res = alpha_hNaF(V)    % NaF inactivation [Cunningham SI 2004].
  hinf = 1.0 ./ (1.0 + exp( (V+58.3)/6.7 ));
  tauh = 0.225 + 1.125 ./ (1.0+exp((V+37.0)/15.0));
  res = hinf./tauh;
end

function res = beta_hNaF(V)    % NaF inctivation  [Cunningham SI 2004].
  hinf = 1.0 ./ (1.0 + exp( (V+58.3)/6.7 ));
  tauh = 0.225 + 1.125 ./ (1.0+exp((V+37.0)/15.0));
  res = (1.0 - hinf)./tauh;
end

function res = alpha_mKDR(V)  % KDR activation [Cunningham SI 2004].
  minf = 1.0 ./ (1.0 + exp( (-V-27.0)/11.5 ));
  taum = 0.25 + 4.35 * exp(-abs(-V-10)/10);
  res = minf./taum;
end

function res = beta_mKDR(V)  % KDR activation [Cunningham SI 2004].
  minf = 1.0 ./ (1.0 + exp( (-V-27.0)/11.5 ));
  taum = 0.25 + 4.35 * exp(-abs(-V-10)/10);
  res = (1.0 - minf)./taum;
end


function aK = alpha_KV(V)
a = 0.0189324;
b = -4.18371;
c = 6.42606;
aK = a*(-((V)+b)) ./ (exp(-((V)+b)/c)-1);
end

function bK = beta_KV(V)
d = 0.015857;
e = 25.4834;
bK = d*exp(-(V)/e);
end

function alpha = alpha_mKM(V)  % 	  ;M-current forward rate function [Traub, 2003].
  alpha = 0.02./(1.0 + exp((-20 - V)/5.));
end

function beta = beta_mKM(V)  % 	  ;M-current backward rate function [Traub, 2003].
  beta = 0.01.*exp((-43 - V)/18.);
end

function res = rs_alpha_m(V)  %    ;NaF inactivation H.  [Traub, 2003]
  hinf = 1.0 ./ (1.0 + exp((-V-34.5)/10.0));
  tauh = 0.0225 + 0.1425*exp(-abs(V+26.5)/10);
  res = hinf./tauh;
end

function res = rs_beta_m(V)  %    ;NaF inactivation H.  [Traub, 2003]
  hinf = 1.0 ./ (1.0 + exp((-V-34.5)/10.0));
  tauh = 0.0225 + 0.1425*exp(-abs(V+26.5)/10);
  res = (1.0 - hinf)./tauh;
end

function res = rs_alpha_h(V)  %    ;NaF inactivation H.  [Traub, 2003]
  hinf = 1.0 ./ (1.0 + exp( (V+59.4)/10.7 ));
  tauh = 0.15 + 1.15 ./ (1.0 + exp( (V+33.5)/15.0 ));
  res = hinf./tauh;
end

function res = rs_beta_h(V)  %    ;NaF inactivation H.  [Traub, 2003]
  hinf = 1.0 ./ (1.0 + exp( (V+59.4)/10.7 ));
  tauh = 0.15 + 1.15 ./ (1.0 + exp( (V+33.5)/15.0 ));
  res = (1.0 - hinf)./tauh;
end

function res = rs_alpha_n(V)  %    ;KDR activation M.  [Traub, 2003]
  minf = 1.0 ./ (1.0 + exp( (-V-29.5)/10.0 ));
  taum = 0.25 + 4.35.*exp(-abs(V+10.0)/10.0);
  res = minf./taum;
end

function res = rs_beta_n(V)  %    ;KDR activation M.  [Traub, 2003]
  minf = 1.0 ./ (1.0 + exp( (-V-29.5)/10.0 ));
  taum = 0.25 + 4.35.*exp(-abs(V+10.0)/10.0);
  res = (1.0 - minf)./taum;
end

function res = rs_alpha_mNaP(V)  %    ;NaP activation.  [Traub, 2003]
  hinf = 1.0 ./ (1.0 + exp((-V-48)/10.0));
  tauh = 0.0225 + 0.1425*exp(-abs(V+40)/10);
  res = hinf./tauh;
end

function res = rs_beta_mNaP(V)   %    ;NaP activation.  [Traub, 2003]
  hinf = 1.0 ./ (1.0 + exp((-V-48)/10.0));
  tauh = 0.0225 + 0.1425*exp(-abs(V+40)/10);
  res = (1.0 - hinf)./tauh;
end

function res = rs_alpha_nKNa(V)        % [Brown J Physiol 2008].
    ka = 0.01;
    eta = -0.0044;
    res = ka*exp(eta*V);
end

function res = rs_beta_nKNa(V)        % [Brown J Physiol 2008].
    kb = 0.0164;
    eta = -0.0206;
    res = kb*exp(eta*V);
end