%%--- UNCONSIDERING PHYSICAL SPACE. No saturation model ---%%
% Carmen Ortega Sabater 
% Latest revision: 2020/10/15
%%-------------------------------------------------------%%


%% Create a directory to save plots and data 

clc
close all
clear all 

run('variables.m')

% Loop 
profile on 

for m = 1:replicates
    for k = 1:simulation_steps

            tic    

              % Unconsidering saturation
                totalN = sum( pheno_space( :, k) );

                % Newborn        
    %           % Adjust prep in line with K 
                prep_modified = (1-(totalN/K)).*prep;
                pdeath_modified = (totalN/K).*pdeath;
                
                %Newborn
                born = binormal( pheno_space(:,k), prep_modified);
                dead = binormal( pheno_space(:,k), pdeath_modified );

                % Update population    
                    pheno_space(:,k) = pheno_space(:,k) + born - dead;

                % Switching phenotypes (only newborn may switch phenotypes)
                    switching = binormal( born, pswitch);                
                    sw_forward = binormal( switching, p_direction);
                    sw_back = switching - sw_forward;

                 for i = 1:pheno
                    if i == 1 
                        pheno_space(i, k) = pheno_space(i, k) - sw_forward(i);
                        pheno_space(i+1, k) = pheno_space( i+1, k) + sw_forward(i);
                    else
                        if  i == size( pheno_space, 1 )
                            pheno_space(i, k) = pheno_space(i, k) - sw_back(i);
                            pheno_space(i-1, k) = pheno_space(i-1, k) + sw_back(i);
                        else
                            pheno_space(i, k) = pheno_space(i, k) - switching(i);    
                            pheno_space(i+1, k) = pheno_space(i+1, k) + sw_forward(i);
                            pheno_space(i-1, k) = pheno_space(i-1, k) + sw_back(i);
                        end  
                    end
                 end

              % Fill next generation with previous generation data   
              pheno_space(:, k+1) = pheno_space(:, k);  


           % Tracking the system 
              history_N(k) = sum(pheno_space(:,k));
    %         if mod( k, savingStep ) == 0    
    %         history_totpop(k/savingStep, k) = sum( pheno_space(:, k));
    %         end

             % Finding average rho 

              % Find total population
                totalN = sum( pheno_space( :, k) );

              % Find N for each phenotype
                totalPheno = pheno_space( :, k);
                rho_average = sum(totalPheno.*prep)/totalN;
                meanRho(k) = rho_average;

              % Find total population
                totpop = sum( pheno_space(:,k));

              % Find frequency for each phenotype
                totalPheno = pheno_space( :, k);
                phenoFreq = totalPheno/totalN;
                history_phenoFreq(k, :) = phenoFreq;

            toc

                % Display time for each iteration
                disp(['Iteration: ' num2str(k)])    

    end
    
    rep_totalN( :, m ) = history_N;
    rep_meanRho( :, m ) = meanRho; 
    rep_phenoFreq( :, :, m ) = history_phenoFreq; 
end 

% Make needed transformations
rep_logtotalN = log( rep_totalN );

profile off
profile viewer 