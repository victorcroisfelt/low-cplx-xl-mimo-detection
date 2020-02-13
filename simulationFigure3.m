% Matlab script used to generate Figure 3 of the article:
%
% "Low-Complexity Distributed XL-MIMO for Multiuser Detection"
%
% @author Victor Croisfelt Rodrigues
% @author Abolfazl Amiri
% @author Taufik Abrao
% @author Elisabeth de Carvalho
% @author Petar Popovski
%

%Mobile communication setup
Mrange = [32 64 128];
Krange = {[16 32],[16 32 64],[16 32 64 128]};

%Number of iterations
numIterRange = {...
    [79 265; 71 270; 71 267]...
    [282 976 0; 75 281 0; 75 281 0],...
    [4010 0 0 0; 66 0 0 0; 66 0 0 0]
    };

%Prepare to save simulation results
RZF = zeros(length(Mrange),length(Krange));
rKA = zeros(length(Mrange),length(Krange),3);

%Go through all M values
for m = 1:length(Mrange)
    
    %Extracting
    M = Mrange(m);
    Kvec = cell2mat(Krange(m));
    
    %Go through all K values
    for k = 1:length(Kvec)
        
        %Extracting
        K = Kvec(k);
        
        %Compute RZF computational complexity
        RZF(m,k) = (3*K.^2*M/2)+(3.*K.*M/2)+((K^3-K)/3) + K;
        
        %Go through all different update schedules
        for us = 1:3

            %Extracting
            numIter = numIterRange{m}(us,k);
            
            %Check update schedule
            if us == 1 % power
                
                %Compute rKA computational complexity
                rKA(m,k,us) = M*numIter + 2*M*K;
                
            else % uniform and active-antenna
                
                %Compute rKA computational complexity
                rKA(m,k,us) = M*numIter + M;
                
            end
            
        end
        
    end
    
end

%% Plot

figure;
hold on; box on; grid on;

surf(log2([16 32 64 128]),log2(Mrange),RZF,'FaceAlpha',0.5,'FaceColor','black')

surf(log2([16 32 64 128]),log2(Mrange),rKA(:,:,1),'FaceAlpha',0.5,'FaceColor',[0 0.4470 0.7410])
surf(log2([16 32 64 128]),log2(Mrange),rKA(:,:,2),'FaceAlpha',0.5,'FaceColor',[0.8500 0.3250 0.0980])
surf(log2([16 32 64 128]),log2(Mrange),rKA(:,:,3),'FaceAlpha',0.5,'FaceColor',[0.9290 0.6940 0.1250])

%shading interp

xlabel('$K$')
ylabel('$M^{(s)}$');
zlabel('Receive combining complexity per subarray');

%xaxis: number of users
xticks([4 5 6 7])
xticklabels({'$2^{4}$','$2^{5}$','$2^{6}$','$2^{7}$'})

%yaxis: subarrays' number of antennas
yticks([5 6 7])
yticklabels({'$2^{5}$','$2^{6}$','$2^{7}$'})

%zaxis: receive combining computational complexity
set(gca,'ZScale','log');

view(135,30)
