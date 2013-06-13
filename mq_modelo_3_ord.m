% Input:
% - Fs              : Sample frequency
% - N_samples       : Número de amostras de cada experimento
% - N_experiments   : Número de experimentos
% - in              : Entrada do sistema (matriz -> linhas são experimentos e colunas são dados)
% - out             : Saída do sistema   (matriz -> linhas são experimentos
% e colunas são dados)
% - use_mode        : Usar moda dos parâmetros para gerar G, senão usa
% esperança amostral
% Output:
% - theta_mean      : Média amostral dos parâmetros obtidos
% - covariance      : Covariância dos parâmetros obtidos
% - theta           : Matriz de resposta por experimento
% - G               : Função de transferência com parâmetros theta_mean

function [theta_mean, covariance, theta, G, Y_est] = mq_modelo_4_ord(Fs, N_samples, N_experiments, in, out, use_mode)
% Classe de modelos Gz = (b0 + b1*z^-1 + b2*z^-2 b3*z^-3 ) / (-a1*z^-1 -a2*z^-2 -a3*z^-3 -a4*z^-4)
%TODO
theta_aux=[];
input = in;
output = out;
NN = N_experiments;
N = N_samples;

for i=1:NN
    %CIs
    y4 = 0; y3 = 0; y2 = 0; y1 = 0; y = 0;
    u3 = 0; u2 = 0; u1 = 0; u0 = 0;   
    
    phi = [];
    Y=[];
    
    for k=1:N
        
        u3 = u2;
        u2 = u1;
        u1 = u0;
        u0 = input(i, k);  
        
        y4 = y3;
        y3 = y2;
        y2 = y1;
        y1 = y;
        
        aux=[y1 y2 y3 y4 u0 u1 u2 u3];
        phi = [phi; aux];

        y = output(i, k);
        Y = [Y;y];
    end

    % Método MQ 
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
a1 = resp(1);a2 = resp(2);a3 = resp(3);a4 = resp(4);
b0 = resp(5);b1 = resp(6);b2 = resp(7);b3 = resp(8);
G = z^-1*(b0 + b1*z^-1 + b2*z^-2 + b3*z^-3 ) / (-a1*z^-1 -a2*z^-2 -a3*z^-3 -a4*z^-4);
Y_est = transpose(phi*theta_mean);

end
