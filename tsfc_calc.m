% ------------------------------------------------------------------------
% Calculate the surface temperature and melt based on the Newton-Raphson
% method and surface energy balance modelling (CLEAN ICE)
% ------------------------------------------------------------------------

%% Initialize 

i = 1;

%% Start the loop

for n = 1:no_gridpointsx
    for l = 1:no_gridpointsy
        
        if mask(n,l) < 1    % If no glacier, SMB is NaN
            
            massbal(i) = NaN;
            
        else

            %% Initialization variables
            
            temp_data_tsfc = temp_glacier(:,n,l);
            temp_tsfc = temp_data_tsfc+273.15;
            prec_tsfc = prec_glacier(:,n,l);
            solarslopingsfc_tsfc = solarslopingsfc(:,n,l);
            solid_prec_glacier_tsfc = solid_prec_glacier(:,n,l);
            lw_in_tsfc = lw_in(:,n,l);
            u_tsfc = u(:,n,l);
            p_tsfc = P(:,n,l);
            tau_tsfc = tau_glacier(:,n,l);
            rh_tsfc = rh_AWS;
            dsnow_tsfc = 0;
            w = 0;
            tsnow_sfc = 0;

            %% Initialization loop

            Td = zeros(1, tmax);
            Ts_past = zeros(1, tmax);
            n_iterations = zeros(1, tmax);
            
            Qnet_tsfc = zeros(1, tmax);
            Lnet_tsfc = zeros(1, tmax);
            Lout_tsfc = zeros(1, tmax);
            SHF_tsfc = zeros(1, tmax);
            Qrain_tsfc = zeros(1, tmax);
            LHF_tsfc = zeros(1, tmax); 
            F_tsfc = zeros(1, tmax);
            
            dQnet_tsfc = zeros(1, tmax);
            dLnet_tsfc = zeros(1, tmax);
            dSHF_tsfc = zeros(1, tmax);
            dQrain_tsfc = zeros(1, tmax);
            dLHF_tsfc = zeros(1, tmax);
            dF_tsfc = zeros(1, tmax);
            
            eS_saturated = zeros(1, tmax);
            eZ_saturated = zeros(1, tmax);
            eS = zeros(1, tmax);
            eZ = zeros(1, tmax);
            qS = zeros(1, tmax);
            qZ = zeros(1, tmax);

                %% Surface temperature and surface energy fluxes calculation  
                  
                for i = 1:tmax
        
                   n_iterations(i) = 0; % Set iterations to vary surface temperature
                   Ts_past(i) = 0;
    
                   % Initially assume Ts = temp, for all other time steps assume it's equal to previous Ts
                   if i == 1
                      Td(1,i) = temp_tsfc(i);
                   else
                      Td(1,i) = Td(1,i-1);
                   end
                    
                   if Td(1,i) > 273.15 
                       Td(1,i) = 273.15; % Correct for maximum of melting surface (0 degrees C)
                   end
                 
                   % Compute fluxes normally
        
                   T_diff(i)=(temp_tsfc(i)-Td(1,i));
                   
                   % Albedo
                   
                   if i == 1
                       alb_tsfc(i) = albsnow+(albice-albsnow).*exp(-dsnow_tsfc(i)./dstar_i);
                       solid_prec_glacier_tsfc(i) = 1e-50;
                   elseif i >= 2
                      tsnow_sfc(i) = -1.*(find(solid_prec_glacier_tsfc(1:i)>0,1,'last')-i);
                      albsnow_tsfc(i) = albfirn+(albsnow-albfirn).*exp(-tsnow_sfc(i)./tsnow);
                      alb_tsfc(i) = albsnow_tsfc(i)+(albice-albsnow_tsfc(i)).*exp(-dsnow_tsfc(i)./dstar_i);
                   end
                   
                   % Surface fluxes
        
                   Qnet_tsfc(i) = tau_tsfc(i).*(1-alb_tsfc(i)).*solarslopingsfc_tsfc(i);
                   Lout_tsfc(i) = em_s.*(stf_bltz.*Td(1,i).^4);
                   Lnet_tsfc(i) = lw_in_tsfc(i) - Lout_tsfc(i);
                   SHF_tsfc(i) = rhoa.*cA.*(C_ex_ice).*u_tsfc(i).*T_diff(i);
                   if temp_tsfc(i) > 273.15
                        eS_saturated(i) = 611.*exp(-Lv./R2.*(1./Td(1,i)-1/273.15));
                        eZ_saturated(i) = 611.*exp(-Lv./R2.*(1./(temp_tsfc(i))-1/273.15));
                        eS(i) = eS_saturated(i);
                        eZ(i) = ((rh_tsfc(i)./100)*eZ_saturated(i));
                        qS(i) = (mwratio).*(eS(i)./p_tsfc(i));
                        qZ(i) = (mwratio).*(eZ(i)./p_tsfc(i));
                        LHF_tsfc(i) = rhoa.*Lv.*(C_ex_ice).*u_tsfc(i).*(qZ(i)-qS(i));
                  else
                        eS_saturated(i) = 0;
                        LHF_tsfc(i) = 0;
                        eZ_saturated(i) = 0;
                        eS(i) = 0;
                        eZ(i) = 0;
                        qS(i) = 0;
                        qZ(i) = 0;
                   end
                   if prec_tsfc(i) > 0 && solid_prec_glacier_tsfc(i) == 0
                       Qrain_tsfc(i) = cW.*rhow.*(prec_tsfc(i)./(sechr)).*(temp_data_tsfc(i));
                   end
                   F_tsfc(i) = Qnet_tsfc(i) + Lnet_tsfc(i) + SHF_tsfc(i) + LHF_tsfc(i) + Qrain_tsfc(i);
        
                   % Calculate derivative of fluxes w.r.t surface temperature
        
                   dQnet_tsfc(i) = 0;
                   dLnet_tsfc(i) = -4.*em_s.*(stf_bltz).*Td(1,i).^3;
                   dSHF_tsfc(i) = -1.*rhoa.*cA.*(C_ex_ice).*u_tsfc(i);
                   if temp_tsfc(i) > 273.15
                        dLHF_tsfc(i) = -1.*rhoa.*Lv.*(C_ex_ice).*u_tsfc(i).*((611.*exp(-Lv./R2.*(1./Td(1,i)-1./273.15)).*(Lv./R2.*Td(1,i).^-2)).*(mwratio/(p_tsfc(i))));
                   else
                        dLHF_tsfc(i) = 0;
                   end
                   if prec_tsfc(i) > 0 && solid_prec_glacier_tsfc(i) == 0
                       dQrain_tsfc(i) = -cW.*rhow.*(prec_tsfc(i)./(sechr));
                   end
                   dF_tsfc(i) = dQnet_tsfc(i) + dLnet_tsfc(i) + dSHF_tsfc(i) + dLHF_tsfc(i);
        
                   % Newton-Raphson method to solve for surface temperature
        
                   while abs(Td(1,i)-Ts_past(i)) > 0.01 && n_iterations < 100
            
                     n_iterations(i) = n_iterations(i) + 1;
                     Ts_past(i) = Td(1,i);
                     Td(1,i) = Ts_past(i) - F_tsfc(i)/dF_tsfc(i); 
                                             
                     if (Td(1,i) - Ts_past(i)) > 1 % Max step size is 1 degree C
                         Td(1,i) = Ts_past(i) + 1;
                     elseif (Td(1,i) - Ts_past(i)) < -1
                         Td(1,i) = Ts_past(i) - 1;
                     end
                     
                     if Td(1,i) > 273.15 
                         Td(1,i) = 273.15; % Correct for maximum of melting surface (0 degrees C)
                     end
                        
                     % Compute fluxes normally
            
                     T_diff(i)=(temp_tsfc(i)-Td(1,i));
                     
                     % Albedo
                     
                     if i == 1
                         alb_tsfc(i) = albsnow+(albice-albsnow).*exp(-dsnow_tsfc(i)./dstar_i);
                         solid_prec_glacier_tsfc(i) = 1e-50;
                     elseif i >= 2
                        tsnow_sfc(i) = -1.*(find(solid_prec_glacier_tsfc(1:i)>0,1,'last')-i);
                        albsnow_tsfc(i) = albfirn+(albsnow-albfirn).*exp(-tsnow_sfc(i)./tsnow);
                        alb_tsfc(i) = albsnow_tsfc(i)+(albice-albsnow_tsfc(i)).*exp(-dsnow_tsfc(i)./dstar_i);
                     end
                     
                     % Surface fluxes
            
                     Qnet_tsfc(i) = tau_tsfc(i).*(1-alb_tsfc(i)).*solarslopingsfc_tsfc(i);
                     Lout_tsfc(i) = em_s.*(stf_bltz.*Td(1,i).^4);
                     Lnet_tsfc(i) = lw_in_tsfc(i) - Lout_tsfc(i);
                     SHF_tsfc(i) = rhoa.*cA.*(C_ex_ice).*u_tsfc(i).*T_diff(i);
                     if temp_tsfc(i) > 273.15
                        eS_saturated(i) = 611.*exp(-Lv./R2.*(1./Td(1,i)-1/273.15));
                        eZ_saturated(i) = 611.*exp(-Lv./R2.*(1./(temp_tsfc(i))-1/273.15));
                        eS(i) = eS_saturated(i);
                        eZ(i) = ((rh_tsfc(i)./100)*eZ_saturated(i));
                        qS(i) = (mwratio).*(eS(i)./p_tsfc(i));
                        qZ(i) = (mwratio).*(eZ(i)./p_tsfc(i));
                        LHF_tsfc(i) = Lv.*rhoa.*(C_ex_ice).*u_tsfc(i).*(qZ(i)-qS(i));
                     else
                        eS_saturated(i) = 0;
                        LHF_tsfc(i) = 0;
                        eZ_saturated(i) = 0;
                        eS(i) = 0;
                        eZ(i) = 0;
                        qS(i) = 0;
                        qZ(i) = 0;
                     end
                     if prec_tsfc(i) > 0 && solid_prec_glacier_tsfc(i) == 0
                        Qrain_tsfc(i) = cW.*rhow.*(prec_tsfc(i)./(sechr)).*(temp_data_tsfc(i));
                     end
                     F_tsfc(i) = Qnet_tsfc(i) + Lnet_tsfc(i) + SHF_tsfc(i) + LHF_tsfc(i) + Qrain_tsfc(i);
        
                     % Calculate derivative of fluxes w.r.t surface temperature
                    
                     dQnet_tsfc(i) = 0;
                     dLnet_tsfc(i) = -4.*em_s.*(stf_bltz).*Td(1,i).^3;
                     dSHF_tsfc(i) = -1.*rhoa.*cA.*(C_ex_ice).*u_tsfc(i);
                     if temp_tsfc(i) > 273.15
                        dLHF_tsfc(i) = -1.*rhoa.*Lv.*(C_ex_ice).*u_tsfc(i).*((611.*exp(-Lv./R2.*(1./Td(1,i)-1./273.15)).*(Lv./R2.*Td(1,i).^-2)).*(mwratio/(p_tsfc(i))));
                     else
                        dLHF_tsfc(i) = 0;
                     end
                     if prec_tsfc(i) > 0 && solid_prec_glacier_tsfc(i) == 0
                        dQrain_tsfc(i) = -cW.*rhow.*(prec_tsfc(i)./(sechr));
                     end
                     dF_tsfc(i) = dQnet_tsfc(i) + dLnet_tsfc(i) + dSHF_tsfc(i) + dLHF_tsfc(i);
          
                     % Set maximum iterations to 100
                     
                     if n_iterations == 100 
                         Td(1,i) = (Td(1,i) + Ts_past(i)) / 2;
                     end

                  end
                                         
                 %% Proceed to surface mass balance calculation
                                  
                 % Calculate snow melt
                 
                 if dsnow_tsfc(i) > 0
                        meltsnow(i) = min(dsnow_tsfc(i),(max(0,(sechr.*F_tsfc(i))./(rhow.*lm))));
                        Wsnow(i) = max(0,(w(i)-(ns.*dsnow_tsfc(i))));
                        w(i+1) = w(i) + meltsnow(i) - Wsnow(i);
                        meltice(i) = 0;
                        runoff(i) = Wsnow(i);
            
                 % In the case of glacier ice surface
            
                 elseif dsnow_tsfc(i) == 0
                       meltice(i) = (max(0,(sechr.*F_tsfc(i))./(rhow.*lm)));
                       Wsnow(i) = 0;
                       w(i+1) = 0;
                       meltsnow(i) = 0;
                       runoff(i) = meltice(i);
                 end
        
                 % Update snow depth
        
                 dsnow_tsfc(i+1) = dsnow_tsfc(i) + solid_prec_glacier_tsfc(i) - meltsnow(i);
        
                 % Avoid negative snow depth
        
                 if dsnow_tsfc(i+1) < 0
                      dsnow_tsfc(i+1) = 0;
                 end
        
                 % Avoid that more water is retained then the retention capacity
        
                 if w(i) > ns.*dsnow_tsfc(i)
                      w(i) = ns.*dsnow_tsfc(i);
                 end
        
                 % Calculate final mass balance profile
        
                 massbal(i+1) = massbal(i) + dt*(solid_prec_glacier_tsfc(i)-runoff(i));
                                 
                end
                 
        end
    end
end

clearvars a_C* b_C* c_C* d_C* A_C* S_C* C j SEB T_diff Qnet_tsfc Lnet_tsfc ...
    SHF_tsfc LHF_tsfc Td dQnet_tsfc dLnet_tsfc dSHF_tsfc dLHF_tsfc ...
    dF_tsfc n_iterations Ts_past eS_saturated eZ_saturated dsnow_tsfc ...
    eS eZ temp_data_tsfc temp_tsfc prec_tsfc solarslopingsfc_tsfc ...
    lw_in_tsfc u_tsfc p_tsfc tau_tsfc alb_tsfc rh_tsfc alb_tsfc tau_tsfc ...
    solid_prec_glacier_tsfcr meltsnow Wsnow w meltice w
