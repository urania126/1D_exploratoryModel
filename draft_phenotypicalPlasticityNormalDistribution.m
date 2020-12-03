% Vector normal de M valores, siendo M el número de fenotipos potenciales
% Crea un vector de 1000 valores aleatorios que se extraen de una
% distribución normal con un promedio de 500 y una desviación de 5
a = 5;
b = 500;
y = a.*randn(1000,1) + b;
sum(y)

pswitch = 0.4 
pstay = 1-pswitch; 

% Fixed distribution: parental phenotype doesn't matter. Intermediate
% phenotypes are more probable
a = randfixedsum( floor(pheno/2), 1, p_switch/2, 0, p_switch/2 )
a2 = sort( a, 'descend' )
a1 = sort( a, 'ascend' )
a_full = [ a1; pstay; a2 ]
plot( a_full ) % ESTO NO VALE!

% Prep depends on the parental phenotype 
% No podemos hacerlo exactamente así

makedist( 'Binomial') 

% binopdf(x, n, p) calcula la función de densidad de probabilidad 
% binomial en cada uno de los valores en el uso del número
% correspondiente de ensayos y la probabilidad de éxito para 
% cada ensayo en xnp

% Funciones: makedist() + pdf()
% binopdf(): https://es.mathworks.com/help/stats/binopdf.html
% Número posible de fenotipos de destino 
% Podemos fijar la anchura de la distribución para permitir
% un número determinado de saltos fenotípicos 

% Normal distribution

x = -25:1:25;
mu = 0;
sigma = 25;
y = normpdf(x,mu,sigma)
figure()
plot( x, y, 'r-x', 'Linewidth', 1.5 )
xlabel( '$Phenotype\; of \; destination (i)$', 'FontSize', 18, 'Interpreter','latex' ); 
ylabel('$P_{j,i}$','FontSize', 18, 'Interpreter','latex');

y(25) = 0 % la probabilidad de moverse en el centro es 0
y = y/sum(y)


% figure()
% hold on
% for i = 1:length( mu )
% y = normpdf(x,mu,sigma)
% plot( x, y )    
% end 
% hold off 

% Binomial distribution
saltos = 0:25; % número de saltos permitidos
               % independiente de la dirección: 0:25
p_stay = 0.2; % Probabilidad de no salto. p_switch = 1-p_stay
newborn = 200; 
y = binopdf( saltos, newborn, p_switch ) 
plot( saltos, y )

%% SUCIO 
% Las probabilidades se van a modificar en función del número 
% de células que hayan nacido??
x_index = 1:50
find( x == x1 )

% Encontrar índices que nos interesan 
in_x = zeros( 1, length(x))
for i = 1:length(x) 
in_x(i) = find( x1 == x(i) )    
end

y1(in_x) - y % La única diferencia es dibujar la curva con más o menos puntos
sum( y/sum(y) )
plot( y/sum(y) )

% Check that everything is OK! 
sum( y)
y - y1(x)

% Moving normal distribution

c = parula(51);
x = -25:1:25;
mu = -25:1:25
sigma = 25;

probMat_sigma25 = zeros( length(x), length(mu) )

figure()
hold on
for i = 1:length( mu )
y = normpdf(x,mu(i),sigma)
y = y/sum(y)
probMat_sigma25(:, i) = y;
plot( x, y, '-x', 'Linewidth', 1.5, 'Color', c(i,:)) 
xlabel( '$Phenotype\; of \; destination (i)$', 'FontSize', 18, 'Interpreter','latex' ); 
ylabel('$P_{j,i}$','FontSize', 18, 'Interpreter','latex');
%pause( 0.025 )
%drawnow
end 
hold off

% narrow distribution: sigma = 5

c = parula(51);
x = -25:1:25;
mu = -25:1:25
sigma = 5;

probMat_sigma5 = zeros( length(x), length(mu) )

figure()
hold on
for i = 1:length( mu )
y = normpdf(x,mu(i),sigma);
y = y/sum(y);
probMat_sigma5(:, i) = y;
plot( x, y, '-x', 'Linewidth', 1.5, 'Color', c(i,:)) 
xlabel( '$Phenotype\; of \; destination (i)$', 'FontSize', 18, 'Interpreter','latex' ); 
ylabel('$P_{j,i}$','FontSize', 18, 'Interpreter','latex');
%pause( 0.025 )
%drawnow
end 
hold off

% Changing width (sigma ranging from 1 to 25) 
x = -25:1:25;
mu = 0
sigma = 1:1:25;
figure()
hold on
for i = 1:length( sigma )
y = normpdf(x,mu,sigma(i))
plot( x, y ) 
xlabel( '$Phenotype\; of \; destination (\rho_{i})$', 'FontSize', 18, 'Interpreter','latex' ); 
ylabel('$P(\rho_{i})$','FontSize', 18, 'Interpreter','latex');
end 
hold off 

%% Sum up 

% Generate probability objects
% Parameters
x = -25:1:25;
mu = -25:1:25
sigma_far = 25;
sigma_close = 5;
mu_fixed = 0;
pheno = 51;

% empty structures
probMat_sigma5 = zeros( length(x), length(mu) );
probMat_sigma25 = zeros( length(x), length(mu) );
probMat_fixed = ones(pheno, pheno);

% Fixed probability (independent from phenotype of origin)
y_fixed = normpdf( x, mu_fixed, sigma_far );
y_fixed = y_fixed/sum(y_fixed);
probMat_fixed = probMat_fixed.*y_fixed';

% Short jumps (but further than to nearest neighbours)
for i = 1:length( mu )
    y = normpdf(x,mu(i),sigma_close);
    y = y/sum(y);
    probMat_sigma5(:, i) = y;
end

% Long jumps
for i = 1:length( mu )
    y = normpdf(x,mu(i),sigma_far);
    y = y/sum(y);
    probMat_sigma25(:, i) = y;
end

% Color plots (imagesc) - Probability matrices

figure()
subplot(1,3,1);
imagesc(probMat_fixed);
colorbar
ylabel( '$Phenotype\; of \; destination (i)$', 'FontSize', 18, 'Interpreter','latex' ); 
xlabel('$Phenotype\; of \; origin (j)$', 'FontSize', 18, 'Interpreter','latex');

subplot(1,3,2);
imagesc(probMat_sigma25)
colorbar 
ylabel( '$Phenotype\; of \; destination (i)$', 'FontSize', 18, 'Interpreter','latex' ); 
xlabel('$Phenotype\; of \; origin (j)$', 'FontSize', 18, 'Interpreter','latex');

subplot(1,3,3);
imagesc(probMat_sigma5);
colorbar
ylabel( '$Phenotype\; of \; destination (i)$', 'FontSize', 18, 'Interpreter','latex' ); 
xlabel('$Phenotype\; of \; origin (j)$', 'FontSize', 18, 'Interpreter','latex');

%% Proof of concept

pheno = 50;
simulation_steps = 100;
Ntimetracking = zeros( simulation_steps, pheno ); 
Prep = ones( pheno, 1)*0.4
pswitch = [0.95 0.05; 0.15 0.85] % las filas de esta matriz tienen que sumar 1
N_0 = ones( pheno, 1)*1000

%Rellenar condiciones iniciales
for i=1:pheno
   Ntimetracking(1, :) =  N_0
end 

% Go! Go! Go!
for t = 1:(simulation_steps-1)
    for i = 1:pheno
        Ntimetracking(t+1,i) = sum(pswitch(i,:).*(1+Prep(i))*Ntimetracking(t, i)) %N_i(t+1)
    end
end 
