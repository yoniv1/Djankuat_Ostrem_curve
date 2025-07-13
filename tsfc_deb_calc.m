% ------------------------------------------------------------------------
% Calculate the surface temperature and melt based on the Newton-Raphson
% and Crank-Nicholson method and surface energy balance modelling
% ------------------------------------------------------------------------
 
% Debris properties
rhor = 2600;                 % Density of debris (kg m-3) from Bozhinskiy et al. 1986
cR = 1260;                   % Specific heat capacity of debris (J kg-1 K-1) from Bozhinskiy et al. 1986
phi_deb = 0.43;              % Debris porosity from Bozhinskiy et al. 1986
albdeb = 0.10;               % Debris albedo from Bozhinskiy et al. 1986
k_r = 2.8;                   % Debris whole rock thermal conductivity (W K-1 m-1) from Bozhinskiy et al. 1986
C_ex_deb = 0.004;            % Exchange coefficient turbulent fluxes debris
hd_crit = 0.03;              % Characteristic snow thickness debris (m w.e.)
phi_low = 0.10;              % Porosity of lowest debris layer
em_d = 0.90;                 % Surface emissivity debris

%% Initialize 

i = 1;

%% Start the loop to determine snow cover evolution

for n = 1:no_gridpointsx
    for l = 1:no_gridpointsy
        
        if mask(n,l) < 2    % If no debris, surface temperature and SMB is NaN (2D glacier mask where debris-covered pixels have value of 2)
            
            yearly_deb_deb(i,n,l) = NaN;
            yearly_runoff_deb(i,n,l) = NaN;
            yearly_massbal_deb_ice(i,n,l) = NaN;
            
        else

            %% Initialization variables
            
            temp_data_deb = temp_glacier(:,n,l);                  % Air temperature
            temp_deb = temp_data_deb+273.15;                      % Transform degrees C to K
            prec_deb = prec_glacier(:,n,l);                       % Precipitation
            solarslopingsfc_deb = solarslopingsfc(:,n,l);         % Solar radiation
            solid_prec_glacier_deb = solid_prec_glacier(:,n,l);   % Solid precipitation
            yearly_melt_clean_ice = yearly_meltice_tsfc(:,n,l);   % Clean ice melt
            lw_in_deb = lw_in(:,n,l);                             
            u_deb = u(:,n,l);                                     % Wind speed clean ice
            u_AWS_deb = u_d(:,n,l);                               % Wind speed debris-covered ice
            p_deb = P(:,n,l);                                     % Air pressure
            tau_deb = tau_glacier(:,n,l);                         % Atmospheric transmissivity
            rh_deb = (max(0,rh_AWS));                             % Relative humidity
            w = 0;
            tsnow_sfc = 0;
            fract_cov = deb_cov(n,l);                             % Fractional debris-covered area 2D mask (ranging 0 to 1)

            %% Initialization loop

            Td = zeros(1, tmax);
            Ts_past = zeros(1, tmax);
            n_iterations = zeros(1, tmax);
            
            Qnet_deb = zeros(1, tmax);
            Lnet_deb = zeros(1, tmax);
            Lout_deb = zeros(1, tmax);
            SHF_deb = zeros(1, tmax);
            LHF_deb = zeros(1, tmax);
            Qrain_sfc_deb = zeros(1, tmax);
            F_deb = zeros(1, tmax);
            
            dQnet_deb = zeros(1, tmax);
            dLnet_deb = zeros(1, tmax);
            dSHF_deb = zeros(1, tmax);
            dLHF_deb = zeros(1, tmax);
            dQrain_sfc_deb = zeros(1, tmax);
            dF_deb = zeros(1, tmax);
            em_d_deb = zeros(1, tmax);
            
            eS_saturated = zeros(1, tmax);
            eZ_saturated = zeros(1, tmax);
            eS = zeros(1, tmax);
            eZ = zeros(1, tmax);
            qS = zeros(1, tmax);
            qZ = zeros(1, tmax);
            tdsnow_deb = zeros(1, tmax);
            dsnow_deb = zeros(1, tmax);
            
            % Compute height of each debris layer with depth
            
            debris_thickness = th_deb(n,l)./100;   % th_deb = 2D debris thickness array
            h = debris_thickness / 10;
            
            % Compute various information needed 
            
            N = debris_thickness/h + 1;
            N_iterations = 0;                   % Set iterations to vary debris thickness
            debris_depth = zeros(N, tmax);      % Debris depth
            phi_deb2 = zeros(N, tmax);          % Debris porosity

            % Initialization numerical schemes
            
            a_Crank = zeros(N,tmax);
            b_Crank = zeros(N,tmax);
            c_Crank = zeros(N,tmax);
            d_Crank = zeros(N,tmax);
            A_Crank = zeros(N,tmax);
            S_Crank = zeros(N,tmax);
            
            Qc_deb = zeros(1, tmax);
            dQc_deb = zeros(1, tmax);
            Qm_ice = zeros(1, tmax);
            Q_rain_deb = zeros(N, tmax);
            
            wd_deb = zeros(N,tmax);
            wd_deb_vol = zeros(N,tmax);
            wd_deb_tot_vol = zeros(1,tmax);
            vol_heat_cap_deb = zeros(N,tmax);
            k_eff_deb = zeros(N,tmax);
            C_deb = zeros(N,tmax);
            
            % Note notation in loop: "i-1" refers to the past

                %% Surface temperature and surface energy fluxes calculation  
                  
                for i = 1:tmax % time loop

                    % Porosity with depth
                    
                    debris_depth(2:N-1,i) = linspace(0,debris_thickness,N-2);
                    phi_deb2(2:N-1,i) = max(0,linspace(phi_deb,phi_low,N-2));

                    if i == 1
                        dsnow_deb(i) = 1e-50;    % Artificial infinitesimally small snow depth for first time step
                        solid_prec_glacier_deb(i) = 1e-50;
                    end
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IN THE CASE SNOW IS PRESENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                   % Calculations for a snow surface
                   
                   if dsnow_deb(i) > 0
        
                   n_iterations(i) = 0; % Set iterations to vary snow surface temperature
                   Ts_past(i) = 0;
                       
                   % Initially assume Ts snow = temp air, for all other time steps assume it's equal to previous Ts
                   if i == 1
                      Td(1,i) = temp_deb(i);
                   else
                      Td(1,i) = Td(1,i-1);
                   end
                    
                   if Td(1,i) > 273.15 
                       Td(1,i) = 273.15; % Correct for maximum of melting surface (0 degrees C)
                   end
                 
                   % Compute fluxes normally
        
                   T_diff(i)=(temp_deb(i)-Td(1,i));
                   
                   % Albedo
                                      
                   if i == 1 % First timestep
                       alb_deb(i) = albsnow+(albdeb-albsnow).*exp(-dsnow_deb(i)./hd_crit);
                   elseif i >= 2 % Else
                      tsnow_sfc(i) = -1.*(find(solid_prec_glacier_deb(1:i)>0,1,'last')-i);
                      albsnow_deb(i) = albfirn+(albsnow-albfirn).*exp(-tsnow_sfc(i)./tsnow);
                      alb_deb(i) = albsnow_deb(i)+(albdeb-albsnow_deb(i)).*exp(-dsnow_deb(i)./hd_crit);
                   end
                   
                   % Surface fluxes

                   % Shortwave net
                   Qnet_deb(i) = tau_deb(i).*(1-alb_deb(i)).*solarslopingsfc_deb(i);
                   % Longwave net
                   Lout_deb(i) = em_s.*(stf_bltz.*Td(1,i).^4);
                   Lnet_deb(i) = lw_in_deb(i) - Lout_deb(i);
                   % Sensible heat flux
                   SHF_deb(i) = rhoa.*cA.*(C_ex_ice).*u_deb(i).*T_diff(i);
                   % Latent heat flux
                   if temp_deb(i) > 273.15
                        eS_saturated(i) = 611.*exp(-Lv./R2.*(1./Td(1,i)-1/273.15));
                        eZ_saturated(i) = 611.*exp(-Lv./R2.*(1./(temp_deb(i))-1/273.15));
                        eS(i) = eS_saturated(i);
                        eZ(i) = ((rh_deb(i)./100)*eZ_saturated(i));
                        qS(i) = (mwratio).*(eS(i)./p_deb(i));
                        qZ(i) = (mwratio).*(eZ(i)./p_deb(i));
                        LHF_deb(i) = rhoa.*Lv.*(C_ex_ice).*u_deb(i).*(qZ(i)-qS(i));
                  else
                        eS_saturated(i) = 0;
                        LHF_deb(i) = 0;
                        eZ_saturated(i) = 0;
                        eS(i) = 0;
                        eZ(i) = 0;
                        qS(i) = 0;
                        qZ(i) = 0;
                   end
                   % Rain heat flux
                   if prec_deb(i) > 0 && solid_prec_glacier_deb(i) == 0;% && dsnow_deb(i) == 0
                       Qrain_sfc_deb(i) = cW.*rhow.*(prec_deb(i)./(sechr)).*(temp_data_deb(i));
                   end
                   % Total energy flux
                   F_deb(i) = Qnet_deb(i) + Lnet_deb(i) + SHF_deb(i) + LHF_deb(i) + Qrain_sfc_deb(i);
        
                   % Calculate derivative of fluxes w.r.t surface temperature
                   
                   % Shortwave net
                   dQnet_deb(i) = 0;
                   % Longwave net
                   dLnet_deb(i) = -4.*em_s.*(stf_bltz).*Td(1,i).^3;
                   % Sensible heat
                   dSHF_deb(i) = -1.*rhoa.*cA.*(C_ex_ice).*u_deb(i);
                   % Latent net
                   if temp_deb(i) > 273.15
                        dLHF_deb(i) = -1.*rhoa.*Lv.*(C_ex_ice).*u_deb(i).*((611.*exp(-Lv./R2.*(1./Td(1,i)-1./273.15)).*(Lv./R2.*Td(1,i).^-2)).*(mwratio/(p_deb(i))));
                   else
                        dLHF_deb(i) = 0;
                   end
                   % Rain heat
                   if prec_deb(i) > 0 && solid_prec_glacier_deb(i) == 0;% && dsnow_deb(i) == 0
                       dQrain_sfc_deb(i) = -cW.*rhow.*(prec_deb(i)./(sechr));
                   end
                   % Total energy flux
                   dF_deb(i) = dQnet_deb(i) + dLnet_deb(i) + dSHF_deb(i) + dLHF_deb(i) + dQrain_sfc_deb(i);
        
                   % Newton-Raphson method to solve for snow surface temperature
        
                   while abs(Td(1,i)-Ts_past(i)) > 0.01 & n_iterations < 100
            
                     n_iterations(i) = n_iterations(i) + 1;
                     Ts_past(i) = Td(1,i);
                     Td(1,i) = Ts_past(i) - F_deb(i)/dF_deb(i); 
                                             
                     if (Td(1,i) - Ts_past(i)) > 1 % Max step size is 1 degree C
                         Td(1,i) = Ts_past(i) + 1;
                     elseif (Td(1,i) - Ts_past(i)) < -1
                         Td(1,i) = Ts_past(i) - 1;
                     end
                     
                     if Td(1,i) > 273.15 
                         Td(1,i) = 273.15; % Correct for maximum of melting surface (0 degrees C)
                     end
                        
                     % Compute fluxes normally
            
                     T_diff(i)=(temp_deb(i)-Td(1,i));
                     
                     % Albedo
                     
                     if i == 1 % First timestep
                         alb_deb(i) = albsnow+(albdeb-albsnow).*exp(-dsnow_deb(i)./hd_crit);
                     elseif i >= 2 % Else
                        tsnow_sfc(i) = -1.*(find(solid_prec_glacier_deb(1:i)>0,1,'last')-i);
                        albsnow_deb(i) = albfirn+(albsnow-albfirn).*exp(-tsnow_sfc(i)./tsnow);
                        alb_deb(i) = albsnow_deb(i)+(albdeb-albsnow_deb(i)).*exp(-dsnow_deb(i)./hd_crit);
                     end
                     
                     % Surface fluxes

                     % Shortwave net
                     Qnet_deb(i) = tau_deb(i).*(1-alb_deb(i)).*solarslopingsfc_deb(i);
                     % Longwave net
                     Lout_deb(i) = em_s.*(stf_bltz.*Td(1,i).^4);
                     Lnet_deb(i) = lw_in_deb(i) - Lout_deb(i);
                     % Sensible heat
                     SHF_deb(i) = rhoa.*cA.*(C_ex_ice).*u_deb(i).*T_diff(i);
                     % Latent heat
                     if temp_deb(i) > 273.15
                        eS_saturated(i) = 611.*exp(-Lv./R2.*(1./Td(1,i)-1/273.15));
                        eZ_saturated(i) = 611.*exp(-Lv./R2.*(1./(temp_deb(i))-1/273.15));
                        eS(i) = eS_saturated(i);
                        eZ(i) = ((rh_deb(i)./100)*eZ_saturated(i));
                        qS(i) = (mwratio).*(eS(i)./p_deb(i));
                        qZ(i) = (mwratio).*(eZ(i)./p_deb(i));
                        LHF_deb(i) = Lv.*rhoa.*(C_ex_ice).*u_deb(i).*(qZ(i)-qS(i));
                     else
                        eS_saturated(i) = 0;
                        LHF_deb(i) = 0;
                        eZ_saturated(i) = 0;
                        eS(i) = 0;
                        eZ(i) = 0;
                        qS(i) = 0;
                        qZ(i) = 0;
                     end
                     % Rain heat
                     if prec_deb(i) > 0 && solid_prec_glacier_deb(i) == 0;% && dsnow_deb(i) == 0
                         Qrain_sfc_deb(i) = cW.*rhow.*(prec_deb(i)./(sechr)).*(temp_data_deb(i));
                     end
                     % Total energy flux
                     F_deb(i) = Qnet_deb(i) + Lnet_deb(i) + SHF_deb(i) + LHF_deb(i) + Qrain_sfc_deb(i);
        
                     % Calculate derivative of fluxes w.r.t surface temperature

                     % Shortwave net
                     dQnet_deb(i) = 0;
                     % Longwave net
                     dLnet_deb(i) = -4.*em_s.*(stf_bltz).*Td(1,i).^3;
                     % Sensible heat
                     dSHF_deb(i) = -1.*rhoa.*cA.*(C_ex_ice).*u_deb(i);
                     % Latent heat
                     if temp_deb(i) > 273.15
                        dLHF_deb(i) = -1.*rhoa.*Lv.*(C_ex_ice).*u_deb(i).*((611.*exp(-Lv./R2.*(1./Td(1,i)-1./273.15)).*(Lv./R2.*Td(1,i).^-2)).*(mwratio/(p_deb(i))));
                     else
                        dLHF_deb(i) = 0;
                     end
                     % Rain heat
                     if prec_deb(i) > 0 && solid_prec_glacier_deb(i) == 0;% && dsnow_deb(i) == 0
                         dQrain_sfc_deb(i) = -cW.*rhow.*(prec_deb(i)./(sechr));
                     end
                     % Total energy flux
                     dF_deb(i) = dQnet_deb(i) + dLnet_deb(i) + dSHF_deb(i) + dLHF_deb(i) + dQrain_sfc_deb(i);
          
                     % Set maximum iterations to 100
                     
                     if n_iterations == 100 
                         Td(1,i) = (Td(1,i) + Ts_past(i)) / 2;
                     end

                   end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IN THE CASE SNOW IS ABSENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Now we do the debris surface

                   elseif dsnow_deb(i) == 0 
                      
                    n_iterations(i) = 0; % Set iterations to vary surface temperature
                    Ts_past(i) = 0;
                    Td(N,i) = 273.15; % Ice temperature is 0 degrees C
                        
                  % For first time step after complete snow melt, initially assume Ts = temp_deb, for all other time steps assume it's equal to previous Ts
                     if i == tdsnow_deb(i-1)+1
                        Td(1,i) = temp_deb(i);
                     else
                        Td(1,i) = Td(1,i-1);
                     end

                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HEAT DIFFUSION/CONDIUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
                  % Calculate temperature profile in the debris
                  
                     if i == tdsnow_deb(i-1)+1 % For first time step after complete snow melt, assume initial condition of linear temperature profile in the debris
            
                          Td_gradient = (Td(1,i) - Td(N,i))/debris_thickness;
            
                         for j = 2:1:(N-1)
                            Td(j,i) = Td(1,i) - (j*h)*Td_gradient;
                         end
                
                     else % Perform Crank-Nicholson Scheme
            
                           for j = 2:1:(N-1) % Equation A8 in Reid and Brock (2010)
                               
                             vol_heat_cap_deb(j,i) = rhor.*cR.*(1-phi_deb2(j,i))+(rhow.*cW.*(wd_deb_vol(j,i)) + rhoa.*cA.*(1-wd_deb_vol(j,i))).*(phi_deb2(j,i));
                             k_eff_deb(j,i) = k_r.*(1-phi_deb2(j,i))+(k_w.*(wd_deb_vol(j,i)) + k_a.*(1-wd_deb_vol(j,i))).*(phi_deb2(j,i));
                             C_deb(j,i) = k_eff_deb(j,i)*sechr/(2*(vol_heat_cap_deb(j,i))*h^2);
                    
                             a_Crank(j,i) = C_deb(j,i);
                             b_Crank(j,i) = 2*C_deb(j,i)+1;
                             c_Crank(j,i) = C_deb(j,i);
                    
                             if j == 2 % Equation A9 in Reid and Brock (2010)
                                d_Crank(j,i) = C_deb(j,i)*Td(1,i) + C_deb(j,i)*Td(1,i-1) + (1-2*C_deb(j,i))*Td(j,i-1) + C_deb(j,i)*Td(j+1,i-1);
                             elseif j < (N-1)
                                d_Crank(j,i) = C_deb(j,i)*Td(j-1,i-1) + (1-2*C_deb(j,i))*Td(j,i-1) + C_deb(j,i)*Td(j+1,i-1);
                             elseif j == (N-1)
                                d_Crank(j,i) = 2*C_deb(j,i)*Td(N,i) + C_deb(j,i)*Td(N-2,i-1) + (1-2*C_deb(j,i))*Td(N-1,i-1);
                             end
                                        
                             if j == 2 % Equations A10 and A11 in Reid and Brock (2010)
                                A_Crank(j,i) = b_Crank(j,i);
                                S_Crank(j,i) = d_Crank(j,i);
                             else 
                                A_Crank(j,i) = b_Crank(j,i) - a_Crank(j,i)/A_Crank(j-1,i)*c_Crank(j-1,i);
                                S_Crank(j,i) = d_Crank(j,i) + a_Crank(j,i)/A_Crank(j-1,i)*S_Crank(j-1,i);
                             end
                    
                           end
                
                           for j = N-1:-1:2 % Equation A12 in Reid and Brock (2010)
                    
                                if j == (N-1)
                                    Td(j,i) = S_Crank(j,i)/A_Crank(j,i);
                                else
                                    Td(j,i) = 1/A_Crank(j,i)*(S_Crank(j,i)+c_Crank(j,i)*Td(j+1,i));
                                end
                    
                           end
                     end
        
                   % Assume snow-free surface and compute fluxes normally
        
                   T_diff(i)=(temp_deb(i)-Td(1,i));
                 
                   if i == 1 % First time step
                      alb_deb(i) = albsnow+(albdeb-albsnow).*exp(-dsnow_deb(i)./hd_crit);
                      solid_prec_glacier_deb(i) = 1e-50;
                   elseif i >= 2 % Else
                     tsnow_sfc(i) = -1.*(find(solid_prec_glacier_deb(1:i)>0,1,'last')-i);
                     albsnow_deb(i) = albfirn+(albsnow-albfirn).*exp(-tsnow_sfc(i)./tsnow);
                     alb_deb(i) = albsnow_deb(i)+(albdeb-albsnow_deb(i)).*exp(-dsnow_deb(i)./hd_crit);
                   end

                   % Shortwave net 
                   Qnet_deb(i) = tau_deb(i).*(1-alb_deb(i)).*solarslopingsfc_deb(i);
                   % Longwave net 
                   Lout_deb(i) = em_d_deb(i).*(stf_bltz.*Td(1,i).^4);
                   Lnet_deb(i) = lw_in_deb(i) - Lout_deb(i);
                   % Sensible heat flux
                   SHF_deb(i) = rhoa.*cA.*(C_ex_deb).*u_AWS_deb(i).*T_diff(i);
                   % Latent heat flux 
                   if prec_deb(i) > 0 && dsnow_deb(i) == 0 && temp_deb(i) > 273.15+t_tresh 
                         eS_saturated(i) = 611.*exp(-Lv./R2.*(1./Td(1,i)-1/273.15));
                         eZ_saturated(i) = 611.*exp(-Lv./R2.*(1./(temp_deb(i))-1/273.15));
                         eS(i) = eS_saturated(i);
                         eZ(i) = ((rh_deb(i)./100)*eZ_saturated(i));
                         qS(i) = (mwratio).*(eS(i)./p_deb(i));
                         qZ(i) = (mwratio).*(eZ(i)./p_deb(i));
                         LHF_deb(i) = Lv.*rhoa.*(C_ex_deb).*u_AWS_deb(i).*(qZ(i)-qS(i));
                    else
                         eS_saturated(i) = 611.*exp(-Lv./R2.*(1./Td(1,i)-1/273.15));
                         eZ_saturated(i) = 611.*exp(-Lv./R2.*(1./(temp_deb(i))-1/273.15));
                         eZ(i) = ((rh_deb(i)./100)*eZ_saturated(i));
                         eS(i) = (eZ(i).*((Td(1,i))./(temp_deb(i))));
                         qS(i) = (mwratio).*(eS(i)./p_deb(i));
                         qZ(i) = (mwratio).*(eZ(i)./p_deb(i));
                         LHF_deb(i) = 0; 
                   end
                   % Rain heat flux
                   if prec_deb(i) > 0 && solid_prec_glacier_deb(i) == 0 && dsnow_deb(i) == 0
                        Qrain_sfc_deb(i) = cW.*rhow.*(prec_deb(i)./(sechr)).*(temp_deb(i)-Td(1,i));
                   end
                   % Conductive heat flux
                   Qc_deb(i) = k_eff_deb(2,i)*(Td(2,i) - Td(1,i))/h;
                   % Total energy flux
                   F_deb(i) = Qnet_deb(i) + Lnet_deb(i) + SHF_deb(i) + LHF_deb(i) + Qrain_sfc_deb(i) + Qc_deb(i);
        
                   % Calculate derivative of fluxes w.r.t surface temperature

                   % Shortwave net 
                   dQnet_deb(i) = 0;
                   % Longwave net 
                   dLnet_deb(i) = -4.*em_d_deb(i).*(stf_bltz).*Td(1,i).^3;
                   % Sensible heat flux 
                   dSHF_deb(i) = -1.*rhoa.*cA.*(C_ex_deb).*u_AWS_deb(i);
                   % Latent heat flux 
                   if prec_deb(i) > 0 && dsnow_deb(i) == 0 && temp_deb(i) > 273.15+t_tresh
                         dLHF_deb(i) = -1.*rhoa.*Lv.*(C_ex_deb).*u_AWS_deb(i).*((611.*exp(-Lv./R2.*(1./Td(1,i)-1./273.15)).*(Lv./R2.*Td(1,i).^-2)).*(mwratio/(p_deb(i))));
                    else
                         dLHF_deb(i) = 0;
                   end
                   % Rain heat flux 
                   if prec_deb(i) > 0 && solid_prec_glacier_deb(i) == 0 && dsnow_deb(i) == 0
                        dQrain_sfc_deb(i) = -cW.*rhow.*(prec_deb(i)./(sechr));
                   end
                   % Conductive heat flux 
                   dQc_deb(i) = -k_eff_deb(2,i)/h;
                   % Total energy flux 
                   dF_deb(i) = dQnet_deb(i) + dLnet_deb(i) + dSHF_deb(i) + dLHF_deb(i) + dQrain_sfc_deb(i) + dQc_deb(i);
        
                 % Newton-Raphson method to solve for surface temperature
        
                    while abs(Td(1,i)-Ts_past(i)) > 0.01 & n_iterations < 100
            
                     n_iterations(i) = n_iterations(i) + 1;
                     Ts_past(i) = Td(1,i);
                     Td(1,i) = Ts_past(i) - F_deb(i)/dF_deb(i); 
                                 
                     if dsnow_deb(i) > 0
                         Td(1,i) = 273.15; % In case of snow cover, surface temperature is 0 degrees C
                     end
            
                     if (Td(1,i) - Ts_past(i)) > 1 % Max step size is 1 degree C
                         Td(1,i) = Ts_past(i) + 1;
                     elseif (Td(1,i) - Ts_past(i)) < -1
                         Td(1,i) = Ts_past(i) - 1;
                     end
                                 
                     % Calculate temperature profile in the debris
            
                     if i == tdsnow_deb(i-1)+1 % For t = 0, which is i = 1, assume initial condition of linear temperature profile in the debris
                         
                          Td(1,i) = temp_deb(i);
                          Td_gradient = (Td(1,i) - Td(N,i))/debris_thickness;
            
                         for j = 2:1:(N-1)
                            Td(j,i) = Td(1,i) - (j*h)*Td_gradient;
                         end
                
                     else % Perform Crank-Nicholson Scheme
            
                         for j = 2:1:(N-1) % Equation A8 in Reid and Brock (2010)
                    
                             vol_heat_cap_deb(j,i) = rhor.*cR.*(1-phi_deb2(j,i))+(rhow.*cW.*(wd_deb_vol(j,i)) + rhoa.*cA.*(1-wd_deb_vol(j,i))).*(phi_deb2(j,i));
                             k_eff_deb(j,i) = k_r.*(1-phi_deb2(j,i))+(k_w.*(wd_deb_vol(j,i)) + k_a.*(1-wd_deb_vol(j,i))).*(phi_deb2(j,i));
                             C_deb(j,i) = k_eff_deb(j,i)*sechr/(2*(vol_heat_cap_deb(j,i))*h^2);
                    
                             a_Crank(j,i) = C_deb(j,i);
                             b_Crank(j,i) = 2*C_deb(j,i)+1;
                             c_Crank(j,i) = C_deb(j,i);
                    
                             if j == 2 % Equation A9 in Reid and Brock (2010)
                                d_Crank(j,i) = C_deb(j,i)*Td(1,i) + C_deb(j,i)*Td(1,i-1) + (1-2*C_deb(j,i))*Td(j,i-1) + C_deb(j,i)*Td(j+1,i-1);
                             elseif j < (N-1)
                                d_Crank(j,i) = C_deb(j,i)*Td(j-1,i-1) + (1-2*C_deb(j,i))*Td(j,i-1) + C_deb(j,i)*Td(j+1,i-1);
                             elseif j == (N-1)
                                d_Crank(j,i) = 2*C_deb(j,i)*Td(N,i) + C_deb(j,i)*Td(N-2,i-1) + (1-2*C_deb(j,i))*Td(N-1,i-1);
                             end
                    
                             if j == 2 % Equations A10 and A11 in Reid and Brock (2010)
                                  A_Crank(j,i) = b_Crank(j,i);
                                  S_Crank(j,i) = d_Crank(j,i);
                             else 
                                  A_Crank(j,i) = b_Crank(j,i) - a_Crank(j,i)/A_Crank(j-1,i)*c_Crank(j-1,i);
                                  S_Crank(j,i) = d_Crank(j,i) + a_Crank(j,i)/A_Crank(j-1,i)*S_Crank(j-1,i);
                             end
                    
                         end
                
                         for j = N-1:-1:2 % Equation A12 in Reid and Brock (2010)
                    
                             if j == (N-1)
                                 Td(j,i) = S_Crank(j,i)/A_Crank(j,i);
                             else
                                 Td(j,i) = 1/A_Crank(j,i)*(S_Crank(j,i)+c_Crank(j,i)*Td(j+1,i));
                             end
                         end
                     end
            
                     % Assume snow-free surface and compute fluxes normally
            
                     T_diff(i)=(temp_deb(i)-Td(1,i));
                     
                     if i == 1 % First timestep
                         alb_deb(i) = albsnow+(albdeb-albsnow).*exp(-dsnow_deb(i)./hd_crit);
                         solid_prec_glacier_deb(i) = 1e-50;
                     elseif i >= 2 % Else
                        tsnow_sfc(i) = -1.*(find(solid_prec_glacier_deb(1:i)>0,1,'last')-i);
                        albsnow_deb(i) = albfirn+(albsnow-albfirn).*exp(-tsnow_sfc(i)./tsnow);
                        alb_deb(i) = albsnow_deb(i)+(albdeb-albsnow_deb(i)).*exp(-dsnow_deb(i)./hd_crit);
                     end

                     % Shortwave net
                     Qnet_deb(i) = tau_deb(i).*(1-alb_deb(i)).*solarslopingsfc_deb(i);
                     % Longwave net
                     Lout_deb(i) = em_d_deb(i).*(stf_bltz.*Td(1,i).^4);
                     Lnet_deb(i) = lw_in_deb(i) - Lout_deb(i);
                     % Sensible heat
                     SHF_deb(i) = rhoa.*cA.*(C_ex_deb).*u_AWS_deb(i).*T_diff(i);
                     % Latent heat
                     if prec_deb(i) > 0 && dsnow_deb(i) == 0 && temp_deb(i) > 273.15+t_tresh
                         eS_saturated(i) = 611.*exp(-Lv./R2.*(1./Td(1,i)-1/273.15));
                         eZ_saturated(i) = 611.*exp(-Lv./R2.*(1./(temp_deb(i))-1/273.15));
                         eS(i) = eS_saturated(i);
                         eZ(i) = ((rh_deb(i)./100)*eZ_saturated(i));
                         qS(i) = (mwratio).*(eS(i)./p_deb(i));
                         qZ(i) = (mwratio).*(eZ(i)./p_deb(i));
                         LHF_deb(i) = Lv.*rhoa.*(C_ex_deb).*u_AWS_deb(i).*(qZ(i)-qS(i));
                     else
                         eS_saturated(i) = 611.*exp(-Lv./R2.*(1./Td(1,i)-1/273.15));
                         eZ_saturated(i) = 611.*exp(-Lv./R2.*(1./(temp_deb(i))-1/273.15));
                         eZ(i) = ((rh_deb(i)./100)*eZ_saturated(i));
                         eS(i) = (eZ(i).*((Td(1,i))./(temp_deb(i))));
                         qS(i) = (mwratio).*(eS(i)./p_deb(i));
                         qZ(i) = (mwratio).*(eZ(i)./p_deb(i));
                         LHF_deb(i) = 0;
                    end
                    % Rain heat
                    if prec_deb(i) > 0 && solid_prec_glacier_deb(i) == 0 && dsnow_deb(i) == 0
                        Qrain_sfc_deb(i) = cW.*rhow.*(prec_deb(i)./(sechr)).*(temp_deb(i)-Td(1,i));
                    end
                    % Conductive heat flux
                    Qc_deb(i) = k_eff_deb(2,i)*(Td(2,i) - Td(1,i))/h;
                    % Total energy flux
                    F_deb(i) = Qnet_deb(i) + Lnet_deb(i) + SHF_deb(i) + LHF_deb(i) + Qrain_sfc_deb(i) + Qc_deb(i);
        
                    % Calculate derivative of fluxes w.r.t surface temperature

                     % Shortwave net
                     dQnet_deb(i) = 0;
                     % Longwave net
                     dLnet_deb(i) = -4.*em_d_deb(i).*(stf_bltz).*Td(1,i).^3;
                     % Sensible heat
                     dSHF_deb(i) = -1.*rhoa.*cA.*(C_ex_deb).*u_AWS_deb(i);
                     % Latent heat
                     if prec_deb(i) > 0 && dsnow_deb(i) == 0 && temp_deb(i) > 273.15+t_tresh 
                         dLHF_deb(i) = -1.*rhoa.*Lv.*(C_ex_deb).*u_AWS_deb(i).*((611.*exp(-Lv./R2.*(1./Td(1,i)-1./273.15)).*(Lv./R2.*Td(1,i).^-2)).*(mwratio/(p_deb(i))));
                      else
                         dLHF_deb(i) = 0;
                     end
                     % Rain flux
                     if prec_deb(i) > 0 && solid_prec_glacier_deb(i) == 0 && dsnow_deb(i) == 0
                        dQrain_sfc_deb(i) = -cW.*rhow.*(prec_deb(i)./(sechr));
                     end
                     % Conductive heat flux
                     dQc_deb(i) = -k_eff_deb(2,i)/h;
                     % Total energy flux
                     dF_deb(i) = dQnet_deb(i) + dLnet_deb(i) + dSHF_deb(i) + dLHF_deb(i) + dQrain_sfc_deb(i) + dQc_deb(i);
          
                     % Set maximum iterations to 100
                     
                     if n_iterations == 100 
                         Td(1,i) = (Td(1,i) + Ts_past(i)) / 2;
                     end

                    end
                
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUB-DEBRIS ENERGY FOR MELTING %%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                    % Calculate energy for melting
                    Qm_ice(i) = max(0,(k_eff_deb(N-1,i)*(Td(N-1,i) - Td(N,i))/h));
        
                    if Qm_ice(i) < 0
                        Qm_ice(i) = 0;
                    end
                       
                   end
                                                                                                                                                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE THE VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                 % Save values for the energy balance components

                 yearly_tempsfc_deb(i,n,l) = Td(1,i);
                 yearly_rh(i,n,l) = rh_deb(i);
                 yearly_ta_deb(i,n,l) = temp_deb(i);
                 yearly_qnet_deb(i,n,l) = Qnet_deb(i);
                 yearly_qin_deb(i,n,l) = tau_deb(i).*solarslopingsfc_deb(i);
                 yearly_albedo_deb(i,n,l) = alb_deb(i);
                 yearly_dsnow_deb(i,n,l) = dsnow_deb(i);
                 yearly_lnet_deb(i,n,l) = Lnet_deb(i);
                 yearly_lout_deb(i,n,l) = Lout_deb(i);
                 yearly_lin_deb(i,n,l) = lw_in_deb(i);
                 yearly_shf_deb(i,n,l) = SHF_deb(i);
                 yearly_lhf_deb(i,n,l) = LHF_deb(i);
                 yearly_qrain_sfc_deb(i,n,l) = Qrain_sfc_deb(i);
                 yearly_psi_deb(i,n,l) = F_deb(i);
                 yearly_ez_saturated_deb(i,n,l) = eZ_saturated(i);
                 yearly_es_saturated_deb(i,n,l) = eS_saturated(i);
                 yearly_es_deb(i,n,l) = eS(i);
                 yearly_ez_deb(i,n,l) = eZ(i);
                 yearly_qs_deb(i,n,l) = qS(i);
                 yearly_qz_deb(i,n,l) = qZ(i);
                 yearly_qc_deb(i,n,l) = Qc_deb(i);
                 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MASS BALANCE CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  
                 % Calculate snow melt
                 
                 if dsnow_deb(i) > 0
                        meltsnow(i) = min(dsnow_deb(i),(max(0,(sechr.*F_deb(i))./(rhow.*lm))));
                        Wsnow(i) = max(0,(w(i)-(ns.*dsnow_deb(i))));
                        w(i+1) = w(i) + meltsnow(i) - Wsnow(i);
                        meltice(i) = 0;
                        runoff(i) = Wsnow(i);
            
                 % In the case of debris-covered ice surface
            
                 elseif dsnow_deb(i) == 0
                       meltice(i) = (fract_cov.*(max(0,(sechr.*Qm_ice(i))./(rhow.*lm)))) + ((1-fract_cov).*yearly_melt_clean_ice(i));
                       Wsnow(i) = 0;
                       w(i+1) = 0;
                       meltsnow(i) = 0;
                       runoff(i) = meltice(i);
                 end
        
                 % Update snow depth
        
                 dsnow_deb(i+1) = dsnow_deb(i) + solid_prec_glacier_deb(i) - meltsnow(i);   % Snow depth
                 tdsnow_deb(i) = (find(dsnow_deb(1:i)>0,1,'last'));                         % Last timestep when snow depth > 0
                 tdsnow_deb2(i) = -1.*(find(dsnow_deb(1:i)>0,1,'last')-i);                  % Time passed since last timestep with snow depth > 0 
        
                 % Avoid negative snow depth
        
                 if dsnow_deb(i+1) < 0
                      dsnow_deb(i+1) = 0;
                 end
        
                 % Avoid that more water is retained then the retention capacity
        
                 if w(i) > ns.*dsnow_deb(i)
                      w(i) = ns.*dsnow_deb(i);
                 end
        
                 % Calculate final mass balance profile
        
                 massbal(i+1) = massbal(i) + dt*(solid_prec_glacier_deb(i)-runoff(i));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE THE VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                 yearly_runoff_deb(i,n,l) = runoff(i);
                 yearly_massbal_deb_ice(i,n,l) = massbal(i);
                 yearly_meltsnow_deb(i,n,l) = meltsnow(i);
                 yearly_Wsnow_deb(i,n,l) = Wsnow(i);
                 yearly_meltice_deb(i,n,l) = meltice(i);
                 
                end
                 
        end
    end
end
