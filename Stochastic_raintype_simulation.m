function[M_raintypes_simul]=Stochastic_raintype_simulation(M_Cov_calib, V_raintypes_calib, Cov_simul, nb_simul, simul_type)

%-----------set simulation---------

V_cov_to_use=[1 2 3 4 5 6];
nb_clusters=max(V_raintypes_calib);

%--------normalize covariance matrices--------
%Rq: normalize both covariance matrices with calibration data (i.e. M_Cov_calib)
M_Cov_calib_ini=M_Cov_calib;
[~,sy]=size(M_Cov_calib);
for i=1:sy
    M_Cov_calib(:,i)=(M_Cov_calib(:,i)-nanmean(M_Cov_calib_ini(:,i)))/nanstd(M_Cov_calib_ini(:,i));
end

[~,sy]=size(Cov_simul);
for i=1:sy
    Cov_simul(:,i)=(Cov_simul(:,i)-nanmean(M_Cov_calib_ini(:,i)))/nanstd(M_Cov_calib_ini(:,i));
end

%-----Initialize structure to store results
M_raintypes_simul=NaN(length(Cov_simul(:,1)),nb_simul);

if simul_type==1 %parametric model
    
    %----------------Step 1: Calibration--------------
    CC = bwconncomp(V_raintypes_calib==0);
    for i=1:length(CC.PixelIdxList)
        if length(CC.PixelIdxList{i})>=24*60/10
            nb_day_dry=floor(length(CC.PixelIdxList{i})/(24*60/10));
            nb_ep_intermittency=length(CC.PixelIdxList{i})-nb_day_dry*(24*60/10);
            if nb_ep_intermittency==0
                nb_day_dry=nb_day_dry-1;
                nb_ep_intermittency=length(CC.PixelIdxList{i})-nb_day_dry*(24*60/10);
            end
            nb_ep_intermittency_start=floor(nb_ep_intermittency/2);
            ep_start=CC.PixelIdxList{i}(1);
            V_raintypes_calib(ep_start+nb_ep_intermittency_start:ep_start+nb_ep_intermittency_start+nb_day_dry*(24*60/10))=-1;
        end
        
    end
    
    %estimate length distribution by rain type (incl. dry)
    V_k=zeros(nb_clusters+2,2);
    V_teta=zeros(nb_clusters+2,2);
    V_str_ecdf=struct();
    
    for my_type=-1:nb_clusters
        CC = bwconncomp(V_raintypes_calib==my_type);
        V_length_rain_event=NaN(length(CC.PixelIdxList),1);
        for i=1:length(CC.PixelIdxList)
            V_length_rain_event(i)=length(CC.PixelIdxList{i}); %in time steps
        end
        
        %empirical pdf
        [f,x] = ecdf(V_length_rain_event);
        V_str_ecdf(my_type+2).x=x;
        V_str_ecdf(my_type+2).f=f;
        
        %parametric pdf (gamma distribution) -> parameter inference
        if my_type>=0
            fun = @(param)-1*((param(1)-1)*sum(log(V_length_rain_event)) - length(V_length_rain_event)*log(gamma(param(1))) - length(V_length_rain_event)*param(1)*log(param(2)) - 1/param(2)*sum(V_length_rain_event));
            param0(1) = 4/(skewness(V_length_rain_event)^2);
            param0(2) = mean(V_length_rain_event)/2;
            A=[[-1 0];[0 -1]];
            b=[1e-10;1e-10];
            param = fmincon(fun,param0,A,b);
            V_k(my_type+2,1)=param0(1);
            V_k(my_type+2,2)=param(1);
            V_teta(my_type+2,1)=param0(2);
            V_teta(my_type+2,2)=param(2);
        end
        
    end
    
    %transition matrix
    TM=zeros(nb_clusters+2,nb_clusters+2);
    str_cov=struct();
    
    for i=1:nb_clusters+2
        for j=1:nb_clusters+2
            str_cov(i,j).data=[];
        end
    end
    
    nb_transition=0;
    for i=2:length(V_raintypes_calib)
        if V_raintypes_calib(i)~=V_raintypes_calib(i-1)
            TM(V_raintypes_calib(i-1)+2,V_raintypes_calib(i)+2)=TM(V_raintypes_calib(i-1)+2,V_raintypes_calib(i)+2)+1;
            nb_transition=nb_transition+1;
            str_cov(V_raintypes_calib(i-1)+2,V_raintypes_calib(i)+2).data=[str_cov(V_raintypes_calib(i-1)+2,V_raintypes_calib(i)+2).data; [M_Cov_calib(i,1) M_Cov_calib(i,2) M_Cov_calib(i,3) M_Cov_calib(i,4) M_Cov_calib(i,5) M_Cov_calib(i,6)]];
        end
    end
    
    CC = bwconncomp(V_raintypes_calib==-1);
    for i=1:length(CC.PixelIdxList)
        nb_day_dry=round(length(CC.PixelIdxList{i})/(24*60/10));
        TM(1,1)=TM(1,1)+nb_day_dry;
        nb_transition=nb_transition+nb_day_dry;
        j=1;
        while j<=length(CC.PixelIdxList{i})
            str_cov(1,1).data=[str_cov(1,1).data; [M_Cov_calib(CC.PixelIdxList{i}(j),1) M_Cov_calib(CC.PixelIdxList{i}(j),2) M_Cov_calib(CC.PixelIdxList{i}(j),3) M_Cov_calib(CC.PixelIdxList{i}(j),4) M_Cov_calib(CC.PixelIdxList{i}(j),5) M_Cov_calib(CC.PixelIdxList{i}(j),6)]];
            j=j+24*60/10;
        end
    end
    TM_Baseline=TM/nb_transition;
    
    %---------Step 2: Simulation---------------
    
    for sim=1:nb_simul
        display(strcat(' Sim:',num2str(sim),'/',num2str(nb_simul)))
        %Simulation: non-homogeneous markov renewal process
        current_type=0;
        current_ind=1;
        V_raintypes_simul=NaN(length(Cov_simul(:,1)),1);
        while current_ind<length(Cov_simul(:,1))
            display(strcat('Time step:',num2str(current_ind),'/',num2str(length(Cov_simul(:,1)))))
            if current_type==-1
                my_duration=24*60/10;
            else
                my_duration = ceil(gamrnd(V_k(current_type+2,2),V_teta(current_type+2,2)));
                if current_type==0
                    my_duration=min(my_duration,24*60/10);
                end
            end
            
            V_raintypes_simul(current_ind:current_ind+my_duration)=current_type;
            current_ind=current_ind+my_duration+1;
            
            if current_ind>=length(Cov_simul(:,1))
                break;
            end
            
            V_P=zeros(1,nb_clusters+2);
            
            for i=1:nb_clusters+2
                
                my_gamma=TM_Baseline(current_type+2,i);
                
                if my_gamma==0 || size(str_cov(current_type+2,i).data,1)<2
                    V_P(i)=0;
                else
                    
                    Xt=zeros(1, length(V_cov_to_use));
                    my_mu=zeros(1, length(V_cov_to_use));
                    my_sigma=zeros(length(V_cov_to_use), length(V_cov_to_use));
                    for ind_cov=1:length(V_cov_to_use)
                        Xt(ind_cov)=Cov_simul(current_ind,V_cov_to_use(ind_cov));
                        my_mu(ind_cov)=nanmean(str_cov(current_type+2,i).data(:,V_cov_to_use(ind_cov)));

                        my_sigma(ind_cov,ind_cov)=nanvar(str_cov(current_type+2,i).data(:,V_cov_to_use(ind_cov)));
                    end
                    
                    V_P(i)=my_gamma*exp(-0.5*(Xt-my_mu)*inv(my_sigma)*(Xt-my_mu)');
                    
                end
            end
            
            V_P_norm=V_P/sum(V_P);
            r = mnrnd(1,V_P_norm);
            current_type=find(r==1)-2;
            
            if isempty(current_type)
                current_type=0;
                display('pb current type')
            end
        end
        
        V_raintypes_simul=V_raintypes_simul(1:length(Cov_simul(:,1)));
        V_raintypes_simul(isnan(V_raintypes_simul))=0;
        V_raintypes_simul(V_raintypes_simul==-1)=0;
        M_raintypes_simul(:,sim)=V_raintypes_simul;
    end
    
else %non-parametric model
    
    %---------------Step 1: define training data----------
    Training_data=[V_raintypes_calib, M_Cov_calib(:,V_cov_to_use)];
    Training_data(Training_data(:,1)==-1,1)=0;
    
    %-------------Step 2: run multiple-point simulation (using QS algorithm)
    Simul_grid=[NaN(length(Cov_simul),1), Cov_simul(:,V_cov_to_use)];
    
    
    kernel2=zeros(601,7);
    kernel2(ceil(length(kernel2(:,1))/2),:)=1;
    kernel2(:,1)=1*10;
    
    serverAddress='localhost';
    for sim=1:nb_simul
        data2=g2s('-sa',serverAddress,'-a','qs','-ti',Training_data(:,:),'-di',Simul_grid(:,:),'-dt',[1,0,0,0,0,0,0],'-ki',kernel2(:,:),'-n',64,'-k',1.2,'-j',128);
        V_raintypes_2=data2(:,1);
        
        CC = bwconncomp(V_raintypes_2==0);
        for i=1:length(CC.PixelIdxList)
            if length(CC.PixelIdxList{i})>=24*60/10
                nb_day_dry=floor(length(CC.PixelIdxList{i})/(24*60/10));
                nb_ep_intermittency=length(CC.PixelIdxList{i})-nb_day_dry*(24*60/10);
                if nb_ep_intermittency==0
                    nb_day_dry=nb_day_dry-1;
                    nb_ep_intermittency=length(CC.PixelIdxList{i})-nb_day_dry*(24*60/10);
                end
                nb_ep_intermittency_start=floor(nb_ep_intermittency/2);
                ep_start=CC.PixelIdxList{i}(1);
                V_raintypes_2(ep_start+nb_ep_intermittency_start:ep_start+nb_ep_intermittency_start+nb_day_dry*(24*60/10))=-1;
            end
        end
        
        V_raintypes_simul=V_raintypes_2;
        M_raintypes_simul(:,sim)=V_raintypes_simul;
    end
end

end