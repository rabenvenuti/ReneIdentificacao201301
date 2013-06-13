close all;
clear all;
clc

%%
%Validação

Fs = 1080;
N_samples = Fs;
N_experiments = 10;

theta_real = [0.5; -0.8; 0; 10]
N = Fs;
NN = N_experiments;

%entrada
Ain=1;
input = ones(1,N).*Ain;
% input = rand([1,N]).*Ain;

%ruido branco
%entrada é mult por 10,
Pnoise = -20;
noise = wgn(2*NN,N,Pnoise);

%ruido colorido
% Criando um filtro passa baixas discreto com:
num_filt = [1];
den_filt = [1 .9];
TS = 1;
Hf = tf(num_filt, den_filt, TS);
for i=1:(2*NN)
    colored_noise(i,:) = filter(num_filt, den_filt, noise(i,:));
end

YY = [];
XX = [];
for i=1:NN
    %CI
    y=0;
    y1=0;
    y2=0;
    u0=0;
    u1=0;
    phi = [];
    Y=[];
    
    for k=1:N
        u1=u0;
        u0=input(k); 
        
        y2=y1;
        y1=y;

        aux=[y1 y2 u0 u1];
        y = aux*theta_real;% + colored_noise(i,k);
        
        % Vetor de saídas
        Y = [Y;y];
    end
    
    YY(i, :) = transpose(Y);
    XX(i, :) = input;
end

%%
in = XX;
out = YY;
use_mode = false;

% method = {'mq'; 'vi'; 'emq'};
% model_order = [2;3];

method = {'mq'};
model_order = [4];

plot_descend = false;
plot_hist = false;
plot_in_out_est = true;
plot_fit_each_theta = false;

%%
% Para cada método
vetor_custos = [];
time = linspace(1/Fs, 1, Fs);
% time = 0:(Fs-1);
for k=1:length(method)
    metodo = char(method(k))
    
    % Para cada ordem de modelo
    for m=1:length(model_order)
        ordem_modelo = model_order(m)

        %nomes
        theta_mean  = genvarname(['mean_' metodo '_' num2str(ordem_modelo)]);
        theta_mode  = genvarname(['mode_' metodo '_' num2str(ordem_modelo)]);
        covariance  = genvarname(['cov_' metodo '_' num2str(ordem_modelo)]);
        theta       = genvarname(['theta_' metodo '_' num2str(ordem_modelo)]);
        G           = genvarname(['G_' metodo '_' num2str(ordem_modelo)]);
        Y_est       = genvarname(['Y_est_' metodo '_' num2str(ordem_modelo)]);
        error       = genvarname(['error_' metodo '_' num2str(ordem_modelo)]);
        custo       = genvarname(['custo_' metodo '_' num2str(ordem_modelo)]);
        custo_tot   = genvarname(['total_custo_' metodo '_' num2str(ordem_modelo)]);
        
        % Resolvo por algum método segundo uma classe de modelos
        eval( [ '[' theta_mean ',' covariance ',' theta ',' G '] = ' ...
            metodo '_modelo_' num2str(ordem_modelo) '_ord(Fs, N_samples, N_experiments, in, out, use_mode);'])
        % Avalio a moda
        eval([theta_mode '= mode(' theta ',2);']);
        
        % Avalio o valor de erro e custo comparando com o 1ro experimento
        eval( [Y_est '= transpose(lsim(' G ', in(1,:)));'])
        eval( [error '= (out(1,:) - ' Y_est ');'])
        eval( [custo '= ' error '.*' error ';'])
        eval( [custo_tot '= sum(' custo ');'])
        vetor_custos = [vetor_custos eval(custo_tot)];        

    end
end

%% Aplico métodos e vario os modelos
vetor_custos = [];

% Para cada método
for k=1:length(method)
    metodo = char(method(k))
    
    % Para cada ordem de modelo
    for m=1:length(model_order)
        ordem_modelo = model_order(m)

        %nomes
        theta_mean  = genvarname(['mean_' metodo '_' num2str(ordem_modelo)]);
        theta_mode  = genvarname(['mode_' metodo '_' num2str(ordem_modelo)]);
        covariance  = genvarname(['cov_' metodo '_' num2str(ordem_modelo)]);
        theta       = genvarname(['theta_' metodo '_' num2str(ordem_modelo)]);
        G           = genvarname(['G_' metodo '_' num2str(ordem_modelo)]);
        Y_est       = genvarname(['Y_est_' metodo '_' num2str(ordem_modelo)]);
        error       = genvarname(['error_' metodo '_' num2str(ordem_modelo)]);
        custo       = genvarname(['custo_' metodo '_' num2str(ordem_modelo)]);
        custo_tot   = genvarname(['total_custo_' metodo '_' num2str(ordem_modelo)]);
        
        % Resolvo por algum método segundo uma classe de modelos
        eval( [ '[' theta_mean ',' covariance ',' theta ',' G '] = ' ...
            metodo '_modelo_' num2str(ordem_modelo) '_ord(Fs, N_samples, N_experiments, in, out, use_mode);'])
        % Avalio a moda
        eval([theta_mode '= mode(' theta ',2);']);
        
        % Avalio o valor de erro e custo comparando com o 1ro experimento
        eval( [Y_est '= transpose(lsim(' G ', in(1,:)));'])
        eval( [error '= (out(1,:) - ' Y_est ');'])
        eval( [custo '= ' error '.*' error ';'])
        eval( [custo_tot '= sum(' custo ');'])
        vetor_custos = [vetor_custos eval(custo_tot)];        

    end
end

% Para cada método
for k=1:length(method)
    metodo = char(method(k));
    
    % Para cada ordem de modelo
    for r=1:length(model_order)
        ordem_modelo = model_order(r);

        %nomes
        theta_mean  = genvarname(['mean_' metodo '_' num2str(ordem_modelo)]);
        covariance  = genvarname(['cov_' metodo '_' num2str(ordem_modelo)]);
        theta       = genvarname(['theta_' metodo '_' num2str(ordem_modelo)]);
        G           = genvarname(['G_' metodo '_' num2str(ordem_modelo)]);
        Y_est       = genvarname(['Y_est_' metodo '_' num2str(ordem_modelo)]);
        error       = genvarname(['error_' metodo '_' num2str(ordem_modelo)]);
        custo       = genvarname(['custo_' metodo '_' num2str(ordem_modelo)]);
        custo_tot   = genvarname(['total_custo_' metodo '_' num2str(ordem_modelo)]);
        
        % Mostro o custo total
        eval(custo_tot)

        % Mostro a cov amostral
        eval(covariance)
        
        % Mostro a moda
%         eval(theta_mode)
        
        % Apresentação dos resultados
        a1_des = sort(eval([theta '(1,:)']), 'descend');
        a2_des = sort(eval([theta '(2,:)']), 'descend');
        
        if model_order(r)==2
            b0_des = sort(eval([theta '(3,:)']), 'descend');
            b1_des = sort(eval([theta '(4,:)']), 'descend');
        else
            a3_des = sort(eval([theta '(3,:)']), 'descend');
            b0_des = sort(eval([theta '(4,:)']), 'descend');
            b1_des = sort(eval([theta '(5,:)']), 'descend');
            b2_des = sort(eval([theta '(6,:)']), 'descend');
        end
        
        % gráfico de amostras descendente
        if plot_descend
            figure
            p=0;
            n=1;
            if model_order(r)==2
                m=4;
                subplot(m,n,p+1), plot(0:length(a1_des)-1, a1_des, 'linewidth', 2); grid on; ylabel('a1'); ...
                    title(['Parâmetros em ordem decrescente para ' metodo ' modelo de ordem ' num2str(ordem_modelo)])
                subplot(m,n,p+2), plot(0:length(a2_des)-1, a2_des, 'linewidth', 2); grid on; ylabel('a2');
                subplot(m,n,p+3), plot(0:length(b0_des)-1, b0_des, 'linewidth', 2); grid on; ylabel('b0');
                subplot(m,n,p+4), plot(0:length(b1_des)-1, b1_des, 'linewidth', 2); grid on; ylabel('b1');
            else
                m=6;
                subplot(m,n,p+1), plot(0:length(a1_des)-1, a1_des, 'linewidth', 2); grid on; ylabel('a1'); ...
                    title(['Parâmetros em ordem decrescente para ' metodo ' modelo de ordem ' num2str(ordem_modelo)])
                subplot(m,n,p+2), plot(0:length(a2_des)-1, a2_des, 'linewidth', 2); grid on; ylabel('a2');
                subplot(m,n,p+3), plot(0:length(a3_des)-1, a3_des, 'linewidth', 2); grid on; ylabel('a3');
                subplot(m,n,p+4), plot(0:length(b0_des)-1, b0_des, 'linewidth', 2); grid on; ylabel('b0');
                subplot(m,n,p+5), plot(0:length(b1_des)-1, b1_des, 'linewidth', 2); grid on; ylabel('b1'); 
                subplot(m,n,p+6), plot(0:length(b2_des)-1, b2_des, 'linewidth', 2); grid on; ylabel('b2');
            end
        end

        % gráficos em histogramas
        if plot_hist
            figure
            aux1=transpose(eval(theta));
            p=0;
            hist_k = 15;
            if model_order(r)==2
                m=2;
                n=2;
                subplot(m,n,p+1), hist(aux1(:,1),hist_k); grid on; ylabel('a1'); 
                subplot(m,n,p+2), hist(aux1(:,2),hist_k); grid on; ylabel('a2');
                    title(['Parâmetros em histograma para ' metodo ' modelo de ordem ' num2str(ordem_modelo)])
                subplot(m,n,p+3), hist(aux1(:,3),hist_k); grid on; ylabel('b0');
                subplot(m,n,p+4), hist(aux1(:,4),hist_k); grid on; ylabel('b1');
            else
                m=2;
                n=3;                
                subplot(m,n,p+1), hist(aux1(:,1),hist_k); grid on; ylabel('a1'); 
                subplot(m,n,p+2), hist(aux1(:,2),hist_k); grid on; ylabel('a2');
                    title(['Parâmetros em histograma para ' metodo ' modelo de ordem ' num2str(ordem_modelo)])
                subplot(m,n,p+3), hist(aux1(:,3),hist_k); grid on; ylabel('a3');
                subplot(m,n,p+4), hist(aux1(:,4),hist_k); grid on; ylabel('b0');
                subplot(m,n,p+5), hist(aux1(:,5),hist_k); grid on; ylabel('b1');
                subplot(m,n,p+6), hist(aux1(:,6),hist_k); grid on; ylabel('b2');
            end
            
        end
       
        % Sistema real e simulado
        if plot_in_out_est
            figure
            plot(time(1, :), in(1,:), time(1, :), out(1,:), time(1, :), eval([Y_est '(1,:)']), 'linewidth', 2)
            grid on;
            title(['Simulação para ' metodo ' modelo de ordem ' num2str(ordem_modelo)]);
            legend('Entrada', 'Saída', 'Estimado com theta médio');
        end

        % Fit para cada metodo
        if plot_fit_each_theta
            z = tf('z');
            tmp = eval(theta);
            [m n] = size(tmp)
            for i=1:n
                resp = eval([theta '(:,i);']);
                if model_order(r)==2
                    a1 = resp(1);a2 = resp(2);
                    b0 = resp(3);b1 = resp(4);
                    GG = (b0*z + b1) / (z*z -a1*z -a2);
                else
                    a1 = resp(1); a2 = resp(2); a3 = resp(3);
                    b0 = resp(4); b1 = resp(5); b2 = resp(6);
                    GG = (b0*z^2 + b1*z + b2) / (z^3 -a1*z^2 -a2*z -a3);
                end
                figure
                plot(time(i,:), out(i, :), time(i,:), transpose(lsim(GG, in(i,:), time(i, :)))) 
                title(['Fit para ' metodo ' modelo de ordem ' num2str(ordem_modelo) ' arquivo ' num2str(i)]);
            end
        end
    end
end
