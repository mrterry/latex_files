function dPsi = derivs(t,Psi); 
% derivs computes the state derivatives for an induction machine 
% The machine is in the synchronous reference frame. 
% Eric Benedict, Spring 1996 
 
% Machine parameters... 
global Rs Rqr Rdr Xls Xlr Xm Xms J P We Wb Vs 
 
% rename state vector... 
Phqs=Psi(1); 
Phds=Psi(2); 
Phqr=Psi(3); 
Phdr=Psi(4); 
Wr=Psi(5); 
Theta=Psi(6); 
 
 
% Select Reference Frame... 
% Assume no angle differnce between a-axis and q-axis at t=0 

% Synchronous
W=We; 
ThetaR=W*t; 
  
% Compute Machine Currents... 
Phmq=Xms*Phqs/Xls + Xms*Phqr/Xlr; 
Phmd=Xms*Phds/Xls + Xms*Phdr/Xlr; 
 
Iqs=(Phqs - Phmq)/Xls; 
Ids=(Phds - Phmd)/Xls; 
Iqr=(Phqr - Phmq)/Xlr; 
                                    Idr=(Phdr - Phmd)/Xlr; 
 
 
% Compute Torques... 
Te=3*P*(Phds*Iqs-Phqs*Ids)/(4*Wb);    %  Electrical Torque 
Tl=0;                                 %  Load Torque 
  
 
% Compute Voltages... 
Vqs=Vs*cos(We*t - ThetaR); 
Vds=-1*Vs*sin(We*t - ThetaR); 
Vqr=0; 
Vdr=0; 
  
% Compute new derivative of state vector... 
dPsi(1) = Vqs - Rs*(Phqs - Phmq)/Xls - W*Phds/Wb;          % PHqs/Wb 
dPsi(2) = Vds - Rs*(Phds - Phmd)/Xls + W*Phqs/Wb;          % PHds/Wb 
dPsi(3) = Vqr - Rqr*(Phqr - Phmq)/Xlr - (W-Wr)*Phdr/Wb;    % PHqr/Wb 
dPsi(4) = Vdr - Rdr*(Phdr - Phmd)/Xlr + (W-Wr)*Phqr/Wb;    % PHdr/Wb 
 
dPsi(5) = (Te - Tl)*P/(2*J*Wb);           % Wr/Wb 
 
dPsi(6) = Wr/Wb;                          % Theta/Wb 

dPsi=Wb*dPsi; 
