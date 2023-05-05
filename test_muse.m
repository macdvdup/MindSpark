clear all; close all;

mu = MuseUdp();

tic
cont = true;
vector_tp9 = [];
vector_af7 = [];
vector_af8 = [];
vector_tp10 = [];
vector_raux = [];
while cont
    [data, timestamp, success] = mu.get_eeg_sample();
    if success
        vector_tp9(end+1) = data(1);
        vector_af7(end+1) = data(2);
        vector_af8(end+1) = data(3);
        vector_tp10(end+1) = data(4);
        vector_raux(end+1) = data(5);
    end
    
    %t=linspace(0,length(vector_tp9),length(vector_tp9));
    %plot(t,vector_tp9)
    %drawnow
    if toc > 65
        cont = false;
    end
end

M = [vector_tp9, vector_af7, vector_af8, vector_tp10, vector_raux];
writematrix(M, 'data_andreia_2.csv')