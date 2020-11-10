%                     \\\\-----VARIABLES-----////
% Carmen Ortega Sabater
% Last review: 14/5/2020 

%28/5/2020
% We need slower tumour growth. We lower prep by a factor of 4. Rho_min
% goes to a value of 0.025 and rho_max to 0.0125
% 101 voxels represent 10cm diameter in reality (more than enough)

% 1/6/2020
% Fixed initial conditions 


%-------------------------------------------------------------------------%

% Seed
%rng( 'shuffle' )
%rng( 1 )

% Time
simulation_steps = 730; % 1 sim_step = 1 cell division = 8h? 
tau = 48; % 1 duplication each 48 hours
day = 24; % hours per day
% ndiv = day/tau; % number of divisions per day
% years = 3 % years of simulation


limit = 5e11; % max number of cells allowed

% Skeleton
    % Number of phenotypes
        pheno = 51;
    % Carrying capacity
        K = 5e09;
    % Number of replicates
        replicates = 50; 

% Skills
    % Proliferation rate
    rho_min = 0.0132*10;
    rho_max = 0.0465*10;
    gap = (rho_max-rho_min)/(pheno-1);
    prep = (rho_min:gap:rho_max);
    prep = (day/tau)*(rho_min:gap:rho_max);
    % NOTE: La probabilidad de que se produzca
    % AL MENOS UNA división en nt horas es 1-(1-p)^n ~ np si p<< 1
    prep = prep'; 
    
    % Phenotype switching
    pswitch_value = 0.4; 
    pswitch = pswitch_value*ones(pheno,1);
    % Allowing switches only to the neighbourhood (equally likely)
    p_direction= 0.5*ones(pheno,1);

    % Death rate (natural)
    pdeath = 0.0298;
    pdeath = pdeath*ones( 1, pheno); 
    %pdeath = 0.25*prep;
  

% Considering therapy
    % Number of cells to begin therapy
        NstartTherapy = 1e9;
    % New needed variables 
    % Define treatment_days vector
        ttdaysPerCycle = 5;                          % Dose days
        ncycles = 5;                                 % Number of cycles
        numberTreatmentDays = ttdaysPerCycle*ncycles; % Dose days per cycle
        treatment_days = zeros( 1, numberTreatmentDays );  

    % Therapy Schemes
        % 23 + 5
        treatmentScheme = [1:5 29:33 57:61 85:89 113:117]; 
         
    % Treatment structure and duration
        nCycles = 5; 
        therapyDays = 5; 
        breakDays = 23; 
        therapyLength = nCycles*(therapyDays + breakDays ); 
    
    % Treatment associated death probabilities
        % Cytotoxic therapy
        pdeath_therapy1 = 2*prep;
        % anti-VEGFR like therapy
        fracDeath_therapy2 = 0.5;
        
% Considering spatial migration
    % Pmig fixed
%      pmigindex = round( length(prep)/2 );
%      pmig = prep(pmigindex)*4;
     %pmig = 0.12/3.5; 
     %pmig = ones(pheno, 1) * pmig;
     
     % Go or grow (update) 
     pmig = linspace(0.01,0.06,pheno);
     pmig = sort( pmig, 'descend');
     
     % What happens if we increase its value 10 times?
     %pmig = 10*pmig; 
    % Go and grow
    % pmig = prep;  
    % Go or grow
    %pmig = sort( prep, 'descend' )
    % Grow and random go
    %pmig = (rho_max-rho_min).*rand(pheno,1)+rho_min  
    
    % We will use also pdirection within spatial migration (?)
        
    % Number of voxels (space)
        dimSpace = 101; 
        
    % K per voxel
        K_voxel = 2e5;

% Tracking the system
    % Phenotype index
        phenoIndex = 1:pheno;
    % Time step Index
        stepIndex = 1:simulation_steps;
    % Save information every X time steps    
        savingStep = 1;
    % AverageRho History
        history_meanRho = zeros( simulation_steps, 1 );
        meanRho = zeros( simulation_steps, 1 );
    % Phenotypic frequencies History
        history_phenoFreq = zeros( simulation_steps, pheno );
    % Population history (all time steps)    
        history_N = zeros( simulation_steps, 1 );
    % Total population by phenotypes considering the whole system at each
    % time step
        totalPheno = zeros( pheno, 1);
    % Phenotypic frequencies considering the whole system at each time step
        phenoFreq = zeros( pheno, 1);
        
    % AverageRho History (considering spatial migration)
    history_meanRho_mig = zeros( simulation_steps, dimSpace );
    % Total population history
    history_totpop_mig = zeros( simulation_steps/savingStep, dimSpace );
    history_totpop_mig_time = zeros( simulation_steps, dimSpace );
    % Space occupation history
    history_space = cell( 1, simulation_steps );
    % meanRho in each voxel 
    meanRho_voxel = zeros( 1, dimSpace );
    
% Track (rep)licates
% these don't require physical space: 

    % Total population for each time step %OK
    rep_totalN = zeros( simulation_steps, replicates);

    % Average rho (taking into account all spatial voxels) for each time
    % step %OK
    rep_meanRho = zeros( simulation_steps, replicates);
    
    % Average phenoFreq total space %OK
    rep_phenoFreq = zeros( simulation_steps, pheno, replicates);
    
% these do require physical space:

    % Average rho per voxel and time step
    rep_meanRhoVoxel = zeros( simulation_steps, dimSpace, replicates); 
    
    % Phenofrequencies per voxel and timestep
    rep_phenoFreqVoxel = zeros( pheno, dimSpace, simulation_steps, replicates);

    % history_totpop_mig (replicates)
    rep_history_totpop_mig = zeros( simulation_steps, dimSpace, replicates);  
    
    % history_meanRho_mig (replicates)
    rep_history_meanRho_mig = zeros( simulation_steps, dimSpace, replicates); 

% Fill initial space
    % No spatial migration
        pheno_space = zeros( pheno, simulation_steps );  
        pheno_start = (pheno+1)/2;
        % Starting with 500 cells with intermediate phenotype
        pheno_space( pheno_start, 1 ) = 500;
        
    % Considering spatial migration
      space_start = (dimSpace+1)/2; %spatial voxel in which we are going to start the simulation
      space = zeros( pheno, dimSpace );   
      
% %   %1. Starting with cells in central voxel (migration considered)
      n_init_1 = 500;     
      space(pheno_start, space_start) = n_init_1;  

%     %2. Starting with a few cells of each phenotype 
%     n_init_2 = ones( pheno, 1)*50; 
%     space(:, space_start) = n_init_2;
    
      %3. Starting with cells following normal distribution
%     n_init_3 = abs( normrnd( 500, 490, [51, 1]) );
%     n_init_3 = floor( sort( n_init_3, 'ascend' ) );  
%     n_ascend = n_init_3(1:25);
%     n_center = n_init_3(26);
%     n_descend = sort( n_ascend, 'descend' );
%     start_normal = [n_ascend; n_center; n_descend];
%     space(:, space_start) = start_normal;