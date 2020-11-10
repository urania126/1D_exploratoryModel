clc
close all
clear all 

%%----------------------SAFE V4------------------------%%
% Carmen Ortega Sabater 
% In this version, we allow prep to become negative to 
% avoid making it 0. 
% Last review: 13/10/2020

% Estimated time of simulation: 18s 

%%-----------------------------------------------------%%

% Run variables script
cd '/Users/carmenortega/Dropbox/PhD/Code/2020_1DEvoDyn_Code/1D_PhysicalSpace_Saturation_NoTherapy'

run('variables.m')

% Track profile
profile on

% Check content of a file
% edit test(.m) 


for m = 1:replicates 
% We have to reset the system before starting every replicate    
% We should include initial conditions in a separate script and run that script here    
    space = zeros( pheno, dimSpace );   
      
  % 1. Starting with cells in central voxel (migration considered)
    n_init_1 = 500;     
    space(pheno_start, space_start) = n_init_1;  

    for k = 1:simulation_steps
           
        tic    
        for j = 1:size(space, 2)

            % Calculating average rho 
            
                % Find total population
                totalN = sum( space, 'all');
                % Find N for each phenotype
                totalPheno = sum( space, 2);
                rho_average = (sum(totalPheno.*prep))/totalN;
                meanRho(k) = rho_average;
            
            % Find total population for each voxel (K is specified per voxel)
                totpop = sum( space(:,j));
                
            % Newborn

            % Considering saturation
                % Adjust prep in line with K 
                prep_modified = (1-(totpop/K_voxel))*prep;
                %pdeath_modified = (totpop/K_voxel)*pdeath;
                
                %Newborn
                born = binormal( space(:,j), prep_modified);
                dead = binormal( space(:, j), pdeath );
                
                % Update population
                space(:,j) = space(:,j) + born - dead;
             
            % Switching phenotypes (only newborn may switch phenotypes)
                switching = binormal( born, pswitch);                
                sw_forward = binormal( switching, p_direction);
                sw_back = switching - sw_forward;
            
             for i = 1:size(space, 1)
                if i == 1 
                    space(i, j) = space(i, j) - sw_forward(i);
                    space(i+1, j) = space( i+1, j ) + sw_forward(i);
                else
                    if  i == size( space, 1 )
                        space(i, j) = space(i, j) - sw_back(i);
                        space(i-1, j) = space(i-1, j) + sw_back(i);
                    else
                        space(i, j) = space(i, j) - switching(i);    
                        space(i+1, j) = space(i+1, j) + sw_forward(i);
                        space(i-1, j) = space(i-1, j) + sw_back(i);
                    end  
                end
             end 
            
            
            % Migration (all cells may migrate)
            
            % Find toppop for j voxel 
            totpop = sum(space(:,j));
                % Adjust pmig in line with K
            pmig_modified = (totpop/K_voxel)*pmig;
            migrating = binormal(space(:,j), pmig_modified);
            mig_forward = binormal(migrating, p_direction);
            mig_back = migrating - mig_forward;
            if j == 1 
                space(:, j) = space(:, j) - mig_forward;
                space(:, j+1) = space(:, j+1) + mig_forward;
            else
                if j == size( space, 2 ) 
                    space(:, j) = space(:, j) - mig_back;
                    space(:, j-1) = space(:, j-1) + mig_back;
                else
                    space(:, j) = space(:, j) - migrating;    
                    space(:, j+1) = space(:, j+1) + mig_forward;
                    space(:, j-1) = space(:, j-1) + mig_back;
                end  
            end
            
     
%         % Tracking the system  
%         
         history_space{k} = space;
         history_N(k) = totalN;
            if mod( k, savingStep ) == 0    
            history_totpop_mig(k/savingStep, j) = sum( space(:, j));
            end
%         
%         % Plot evolution phenotypic space
%         
        phenoFreq = totalPheno/totalN;
        history_phenoFreq(k, :) = phenoFreq;
        
%         
        % Plot evolution spatial dimension
        totalN_voxel = sum( space(:,j)); %
        
        meanRho_voxel(j) = (sum(space(:,j).*prep))/sum(space(:,j));
    
        
        % Pheno frequency per voxel
        voxelPhenoFreq = space(:,j)/sum(space(:,j)); 
        
            % Track phenotypic frequencies per voxel and time step    
        rep_phenoFreqVoxel(:,j,k,m) = voxelPhenoFreq; %OK;
        % Dim 1 - Fenotipos 
        % Dim 2 - Voxel 
        % Dim 3 - Tiempo. ¿Es necesario guardarlo para todos los puntos
        % temporales?
        % Dim 4 - Réplicas 
        
        
            % Esto a veces genera NaN porque estamos dividiendo 0/0 en los
            % voxels que aún no han sido invadidos 
                toc
        
        % Display time for each iteration
        disp(['Iteration: ' num2str(k)]) 
        end
   
    history_meanRho_mig( k, : ) = meanRho_voxel; %
    
    % Track phenotypic frequencies for the whole tumor at each time step
    rep_phenoFreq( k, :, m ) =  phenoFreq;
        % Dim 1: simulation_steps
        % Dim 2: pheno
        % Dim 3: replicates
    
    end
 rep_totalN(:, m) = history_N; %OK
 rep_meanRho(:,m) = meanRho; %OK 
 rep_meanRhoVoxel( :, :, m) = history_meanRho_mig; %OK
 rep_history_totpop_mig(:,:,m) = history_totpop_mig;
 rep_history_meanRho_mig(:, :, m) = history_meanRho_mig; 
 
 
 %rep_phenoFreqVoxel(:,j,k,m) = voxelPhenoFreq; %OK; dentro del bucle k 
    % Dim 1 - Fenotipos 
    % Dim 2 - Voxel 
    % Dim 3 - Tiempo. ¿Es necesario guardarlo para todos los puntos
    % temporales?
    % Dim 4 - Réplicas 
    
    % Y si almacenamos estos datos en un solo objeto tridimensional
    % escogiendo los voxels de los extremos y el voxel central?        
end

profile off
profile viewer 

% export the results

% Make needed transformations
rep_logtotalN = log( rep_totalN );

