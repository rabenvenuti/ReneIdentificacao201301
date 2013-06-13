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

function [theta_mean, covariance, theta, G] = vi_modelo_3_ord(Fs, N_samples, N_experiments, in, out, use_mode)
% Classe de modelos Gz = (b0*z^2 + b1*z + b2) / (z^3 -a1*z^2 -a2*z -a3)

theta_aux=[];
input = in;
output = out;
NN = N_experiments;
N = N_samples;

for i=1:NN
    %CI
    y=0;
    y1=0;
    y2=0;
    y3=0;
    u0=0;
    u1=0;
    u2=0;
    
    phi_aux = [];
    Y=[];
    
    for k=1:N
        u2 = u1;
        u1 = u0;
        u0 = input(i, k);   

        y3 = y2;
        y2 = y1;
        y1 = y;
        
        aux=[y1 y2 y3 u0 u1 u2];
        phi_aux = [phi_aux; aux];
                
        y = output(i,k);
        Y = [Y;y];
    end

    % Usando como instrumento um experimento diferente
    if mod(i,2)==1
        %Se �mpar, salvo o instrumento
        Z= phi_aux;
    else
        %Se par, salvo phi e uso o instrumento adquirido 
        %no experimento anterior para c�lculo dos 
        %par�metros
        phi = phi_aux;

        aux_matrix = Z;
        matrix = phi;
        theta_aux = [theta_aux, inv(transpose(aux_matrix)*matrix)*transpose(aux_matrix)*Y];

        phi= [];
        Z=[];
     end
end

% Theta vector
theta = theta_aux;

% Expected value
theta_mean = transpose(mean(transpose(theta)));

% Moda value
theta_mode = mode(theta, 2);

% Covariance
covariance = diag(cov(transpose(theta)));

Ts = 1/Fs;
z = tf('z',Ts);
if use_mode
    resp = theta_mode;
else
    resp = theta_mean;
end
a1 = resp(1); a2 = resp(2); a3 = resp(3);
b0 = resp(4); b1 = resp(5); b2 = resp(6);
G = (b0*z^2 + b1*z + b2) / (z^3 -a1*z^2 -a2*z -a3);

end
