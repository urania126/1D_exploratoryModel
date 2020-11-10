 % Saving the results

% cd .. % go back to previous folder
% mkdir results
cd ./results
mkdir test31
cd ./test31

save( 'rep_totalN.mat', 'rep_totalN' );
save( 'rep_meanRho.mat', 'rep_meanRho' );
save( 'rep_phenoFreq.mat', 'rep_phenoFreq' );
save( 'rep_logtotalN.mat', 'rep_logtotalN' );

save( 'rep_phenoFreqVoxel.mat', 'rep_phenoFreqVoxel' );
save( 'rep_history_meanRho_mig.mat', 'rep_history_meanRho_mig');
save( 'rep_history_totpop_mig.mat', 'rep_history_totpop_mig'); 

save( 'rep_phenoFreqVoxel_1.mat', 'rep_phenoFreqVoxel_1' );
save( 'rep_phenoFreqVoxel_2.mat', 'rep_phenoFreqVoxel_2' );

% Move to proper directory

% load mean and sd files
t3_meanRho = load('mean_meanRho.mat');
t3_meanRho = t3_meanRho.mean_meanRho;
t3_sdmeanRho = load( 'sd_meanRho.mat');
t3_sdmeanRho = t3_sdmeanRho.sd_meanRho;

% Import the data
rep_meanRho = load( 'rep_meanRho.mat' );
rep_meanRho = rep_meanRho.rep_meanRho;

rep_phenoFreq = load( 'rep_phenoFreq.mat' );
rep_phenoFreq = rep_phenoFreq.rep_phenoFreq;

rep_phenoFreqVoxel = load( 'rep_phenoFreqVoxel.mat'); 
rep_phenoFreqVoxel = rep_phenoFreqVoxel.rep_phenoFreqVoxel;

rep_history_meanRho_mig = load( 'rep_history_meanRho_mig.mat' );
rep_history_meanRho_mig = rep_history_meanRho_mig.rep_history_meanRho_mig;

rep_history_totpop_mig = load( 'rep_history_totpop_mig.mat' );
rep_history_totpop_mig = rep_history_totpop_mig.rep_history_totpop_mig;

rep_totalN = load( 'rep_totalN.mat' );
rep_totalN = rep_totalN.rep_totalN;

rep_logtotalN = load( 'rep_logtotalN.mat' );
rep_logtotalN = rep_logtotalN.rep_logtotalN;

% 1. meanRho file
mean_meanRho = mean( rep_meanRho, 2 ); % mean by rows
sd_meanRho = std(rep_meanRho,[],2); %standard deviation by rows
sem_meanRho = sd_meanRho/sqrt(replicates);

% 2. history_meanRho_mig (meanRhoVoxel) 
mean_history_meanRho_mig = mean( rep_history_meanRho_mig, 3 );
sd_history_meanRho_mig = std( rep_history_meanRho_mig,[],3);
sem_history_meanRho_mig = sd_history_meanRho_mig/sqrt(replicates); 

% 3. phenoFreq 
mean_phenoFreq = mean( rep_phenoFreq, 3 );
sd_phenoFreq = std(rep_phenoFreq,[],3);
sem_phenoFreq = sd_phenoFreq/sqrt(replicates);

% 4. phenoFreqVoxel
index_time = floor([1:simulation_steps/4:simulation_steps, simulation_steps ]);
varname = genvarname({'T', 'T', 'T', 'T', 'T'});

    % Generate objects' names
    cutPhenoFreqVoxel = zeros( pheno, dimSpace, length(index_time), replicates);
    objNaming = cell([length(index_time), 1]);

    % i think these lines are useless but maybe they won't in the future:
    % for t = 1:length(index_time)
    %     timeName = sprintf('phenoFreqVoxel%d', t);
    %     objNaming{t} = timeName; 
    % end

    for t = 1:length(index_time)
        cutPhenoFreqVoxel(:,:,t,:) = rep_phenoFreqVoxel(:,:,index_time(t),:);
    end

mean_phenoFreqVoxel = mean( cutPhenoFreqVoxel, 4); 
sd_phenoFreqVoxel = std( cutPhenoFreqVoxel, [], 4); 
sem_phenoFreqVoxel = sd_phenoFreqVoxel/sqrt(replicates);

% 5. History totpop mig 
mean_history_totpop_mig = mean( rep_history_totpop_mig, 3 );
sd_history_totpop_mig = std( rep_history_totpop_mig,[],3);
sem_history_totpop_mig = sd_history_totpop_mig/sqrt(replicates);

% 6. totalN
mean_totalN = mean( rep_totalN, 2 ); % mean by rows
sd_totalN = std( rep_totalN,[],2); %standard deviation by rows
sem_totalN = sd_totalN/sqrt(replicates);

% 7. Log totalN
mean_logtotalN = mean( rep_logtotalN, 2 ); % mean by rows
sd_logtotalN = std( rep_logtotalN,[],2); %standard deviation by rows
sem_logtotalN = sd_logtotalN/sqrt(replicates);

% Save files 

save( 'mean_meanRho.mat', 'mean_meanRho' );
save( 'sd_meanRho.mat', 'sd_meanRho' );
save( 'sem_meanRho.mat', 'sem_meanRho' );

save( 'mean_history_meanRho_mig.mat', 'mean_history_meanRho_mig' );
save( 'sd_history_meanRho_mig.mat', 'sd_history_meanRho_mig' );
save( 'sem_history_meanRho_mig.mat', 'sem_history_meanRho_mig');

save( 'mean_phenoFreq.mat', 'mean_phenoFreq' );
save( 'sd_phenoFreq.mat', 'sd_phenoFreq' );
save( 'sem_phenoFreq.mat', 'sem_phenoFreq');

save( 'mean_phenoFreqVoxel.mat', 'mean_phenoFreqVoxel' );
save( 'sd_phenoFreqVoxel.mat', 'sd_phenoFreqVoxel' );
save( 'sem_phenoFreqVoxel.mat', 'sem_phenoFreqVoxel');

save( 'mean_totalN.mat', 'mean_totalN' );
save( 'sd_totalN.mat', 'sd_totalN' );
save( 'sem_totalN.mat', 'sem_totalN');

save( 'mean_logtotalN.mat', 'mean_logtotalN' );
save( 'sd_logtotalN.mat', 'sd_logtotalN' );
save( 'sem_logtotalN.mat', 'sem_logtotalN');

save( 'mean_history_totpop_mig.mat', 'history_totpop_mig' );
save( 'sd_history_totpop_mig.mat', 'sd_history_totpop_mig' );
save( 'sem_history_totpop_mig.mat', 'sem_history_totpop_mig');