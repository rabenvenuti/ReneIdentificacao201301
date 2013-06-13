% Input:
% - Fs              : Sample frequency
% - N_samples       : N�mero de amostras de cada experimento
% - N_experiments   : N�mero de experimentos
% - in              : Entrada do sistema (matriz -> linhas s�o experimentos e colunas s�o dados)
% - out             : Sa�da do sistema   (matriz -> linhas s�o experimentos
% e colunas s�o dados)
% - use_mode        : Usar moda dos par�metros para gerar G, sen�o usa
% esperan�a amostral
% Output:
% - theta_mean      : M�dia amostral dos par�metros obtidos
% - covariance      : Covari�ncia dos par�metros obtidos
% - theta           : Matriz de resposta por experimento
% - G               : Fun��o de transfer�ncia com par�metros theta_mean

function [theta_mean, covariance, theta, G] = mq_modelo_2_ord(Fs, N_samples, N_experiments, in, out, use_mode)
% Classe de modelos Gz = (b0*z + b1) / (z^2 -a1*z -a2)

theta_aux=[];
input = in;
output = out;
NN = N_experiments;
N = N_samples;

for i=1:NN
    %CIs
    y2 = 0; y1 = 0; y = 0;
    u1 = 0; u0 = 0;   
    
    phi = [];
    Y=[];
    
    for k=1:N
        u1 = u0;
        u0 = input(i, k);  
        
        y2 = y1;
        y1 = y;
        
        aux=[y1 y2 u0 u1];
        phi = [phi; aux];

        y = output(i, k);
        Y = [Y;y];
    end

    % M�todo MQ 
    aux_matrix = phi;
    matrix = phi;
    theta_aux = [theta_aux, inv(transpose(aux_matrix)*matrix)*transpose(aux_matrix)*Y];
end

% Theta vector
theta = theta_aux;

% Expected value
theta_mean = transpose(mean(transpose(theta)));

% Moda value
theta_mode = mode(theta, 2);

% Covariance
covariance = cov(transpose(theta));

Ts = 1/Fs;
z = tf('z',Ts);
if use_mode
    resp = theta_mode;
else
    resp = theta_mean;
end
a1 = resp(1);a2 = resp(2);
b0 = resp(3);b1 = resp(4);
G = (b0*z + b1) / (z*z -a1*z -a2);

end
