function [in, out, time] = data_loader(N, Fs, div_factor, initial_file_s, initial_file_nbr, folder)
% Input:
% - N                   : Number of files
% - Fs                  : Sample frequency
% - div_factor          : Over sample divider
% - initial_file_s      : Mat file pré number name
% - initial_file_nbr    : Mat file initial value
% - folder              : Data location folder
% Output:
% - in                  : Input variable
% - out                 : Output variable
% - time                : Sample time vector

aux = zeros(Fs, 1);
cont = 1;
init = initial_file_nbr;

for i = 1:N
    
    filename = sprintf('%s%2.2d', initial_file_s, (init + i));
    
    load( [folder filename '.mat'] );

    eval( ['File ='  filename ]);
    
    time(i, :) = linspace(min(File.X.Data), max(File.X.Data), max(File.X.Data)*Fs);
        
    [row numCaptures] = size(File.Y);
    
    for num = 1:numCaptures
        
        tmp = File.Y(num).Data;
        aux = [];
        
        for j = 1:length(tmp)
            
            if (mod(j, div_factor) == 0)
               aux = [aux (sum(tmp((j-(div_factor-1)):j)))/div_factor];
            end
            
        end
        
        if num == 1
            in(i, :) = aux;
        else
            out(i, :) = aux;
        end
        
    end
        
end