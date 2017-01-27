function trans_thread_head(tt1,polno,impolno)
% impolno = 1;
% polno=1;
% tt1=1;

% global Tss T nz nk kgrid nb bgrid

load params.mat
% load Pop_Evol.mat
load Surv_Probs.mat
% load Imm_Data.mat
filename = ['imm_polparams_' num2str(impolno) '.mat'];
load(filename)
% polno
% impolno
% size(im_rate)
% size(imm_age)
% size(pgrowth)
% if length(im_rate)<Tss
%     leg_rate = [leg_rate',leg_rate(end).*ones(1,Tss-(length(leg_rate)-1))];
% end
% 
% if length(im_rate)>=Tss
%     leg_rate2 = [leg_rate(tt1:Tss)', leg_rate(Tss).*ones(1,T-(Tss-(tt1-1)))];
% else
%     leg_rate2 = [leg_rate(tt1:Tss), leg_rate(Tss).*ones(1,T-(Tss-(tt1-1)))];
% end
% 
% 
% if length(im_rate)<Tss
%     illegal_rate = [illegal_rate',illegal_rate(end).*ones(1,Tss-(length(illegal_rate)-1))];
% end
% 
% if length(im_rate)>=Tss
%     illeg_rate2 = [illegal_rate(tt1:Tss)', illegal_rate(Tss).*ones(1,T-(Tss-(tt1-1)))];
% else
%     illeg_rate2 = [illegal_rate(tt1:Tss), illegal_rate(Tss).*ones(1,T-(Tss-(tt1-1)))];
% end
% imm_age2 = [imm_age(tt1:Tss), imm_age(end).*ones(1,T-(Tss-(tt1-1)))];
% pgrowth2 = [pgrowth(tt1:Tss), pgrowth(end).*ones(1,T-(Tss-(tt1-1)))];


% MU2 = zeros(ndem,T);
% 
% for demtype = 1:ndem
%     MU2(demtype,:) =  dist_year(1,demtype).*((1/sum(Mu2).*Mu2));
% end


dist_1 = zeros(nk,nz,nb,min(T+1,Tss+1),3,ndem);
dist_r = zeros(nk,nb,min(T+1,Tss+1),3,ndem);

cohort = tt1+(T-1);

% load SSVALS pop_prev
load SSVALS pop_prev
pop_trans = [pop_prev; pop_trans];
% legal_rate = [legal_rate(1); legal_rate];
% illegal_rate = [illegal_rate(1); illegal_rate];


for demtype = 1:ndem
%     demtype_ni = (demtype-1)*1 + (1-(demtype-1))*2;   % gives the other demtype (2 gives 1 and 1 gives 2)
%     demtype = 2;


    filename = ['head' num2str(tt1) '_' num2str(demtype) '_' num2str(polno) '.mat'];
    jobdir = ['job_' num2str(polno) '_' num2str(impolno)];  % gives name
    totfile = fullfile(jobdir,filename);
    load(totfile);


    year = max(1,min(cohort +1 - T,Tss));  % with age = 1
    
%         while beqeps>beqtol
%     im_flow = [pop_trans(year)*(pgr); pop_trans(year)*imm_age(1)*legal_rate(year); pop_trans(year)*imm_age(1)*illegal_rate(year)];   % using period 1 imm rate values for steady state
    im_flow = [pop_trans(year)*(pgr); pop_trans(year)*imm_age(1)*legal_rate(1); pop_trans(year)*imm_age(1)*illegal_rate(1)];   % using period 1 imm rate values for steady state
    loc1 = 1; loc2 = 1;  % initiates k and b at zero.
    

    for i1=1:nz
        for immigrant_type = 1:3
            dist_1(loc1,i1,loc2,1,immigrant_type,demtype) = proddist_age(i1,1,immigrant_type)*im_flow(immigrant_type);
        end
    end
    
    for t1 = 1:min(Tr,Tss-(tt1-1)-1) % Tss-(tt1-1) gives age in final period of transition
        %----------------------------------CHECK--------------------------
        %----------------------------------CHECK--------------------------
        %----------------------------------CHECK--------------------------
        age = t1;
        year = max(1,min(cohort +age - T,Tss))+1;
        %----------------------------------CHECK--------------------------
        %----------------------------------CHECK--------------------------
        %----------------------------------CHECK--------------------------
%         im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(year); pop_trans(year)*imm_age(age)*illegal_rate(year)];   % using period 1 imm rate values for steady state
        im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(1); pop_trans(year)*imm_age(age)*illegal_rate(1)];   % using period 1 imm rate values for steady state
        for j2 = 1:nz
            for i1 = 1:nk
                for i2 = 1:nb
                    point_k = max(kopt(i1,j2,i2,t1),kgrid(1));    % placing floor at lowest gridpoint.
                    loc1 = find(kgrid(1:nk-1)<=point_k,1,'last');
                    w1 = (point_k - kgrid(loc1))/(kgrid(loc1+1)-kgrid(loc1));  % amount allocated to higher gridpoint
                    w1 = min(w1,1);

                    point_b = max(bopt(i1,j2,i2,t1),bgrid(1));    % placing floor at lowest gridpoint.
                    loc2 = find(bgrid(1:nb-1)<=point_b,1,'last');    % lower gridpoint
                    w2 = (point_b - bgrid(loc2))/(bgrid(loc2+1)-bgrid(loc2));  % amount allocated to higher gridpoint
                    w2 = min(w2,1);
                    
                    for j4 = 1:nz
                        
                        if t1<Tr
                            for immigrant_type = 1:3
                                dist_hold = dist_1(i1,j2,i2,t1,immigrant_type,demtype) + eq(i1,1)*eq(i2,1)*proddist_age(j2,t1,immigrant_type)*im_flow(immigrant_type);

                                dist_1(loc1,j4,loc2,t1+1,immigrant_type,demtype) = dist_1(loc1,j4,loc2,t1+1,immigrant_type,demtype) + surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)*dist_hold;
                                dist_1(loc1+1,j4,loc2,t1+1,immigrant_type,demtype) = dist_1(loc1+1,j4,loc2,t1+1,immigrant_type,demtype) + surv(t1)*(1-w2)*(w1)*tr_z(j2,j4)*dist_hold;
                                dist_1(loc1,j4,loc2+1,t1+1,immigrant_type,demtype) = dist_1(loc1,j4,loc2+1,t1+1,immigrant_type,demtype) + surv(t1)*(w2)*(1-w1)*tr_z(j2,j4)*dist_hold;
                                dist_1(loc1+1,j4,loc2+1,t1+1,immigrant_type,demtype) = dist_1(loc1+1,j4,loc2+1,t1+1,immigrant_type,demtype) + surv(t1)*(w2)*(w1)*tr_z(j2,j4)*dist_hold;
                            end
                            
                        
                        elseif t1==Tr
                            for immigrant_type = 1:3
                                dist_hold = dist_1(i1,j2,i2,t1,immigrant_type,demtype) + eq(i1,1)*eq(i2,1)*proddist_age(j2,t1,immigrant_type)*im_flow(immigrant_type);

                                dist_r(loc1,loc2,t1+1,immigrant_type,demtype) = dist_r(loc1,loc2,t1+1,immigrant_type,demtype) + surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)*dist_hold;
                                dist_r(loc1+1,loc2,t1+1,immigrant_type,demtype) = dist_r(loc1+1,loc2,t1+1,immigrant_type,demtype) + surv(t1)*(1-w2)*w1*tr_z(j2,j4)*dist_hold;
                                dist_r(loc1,loc2+1,t1+1,immigrant_type,demtype) = dist_r(loc1,loc2+1,t1+1,immigrant_type,demtype) + surv(t1)*w2*(1-w1)*tr_z(j2,j4)*dist_hold;
                                dist_r(loc1+1,loc2+1,t1+1,immigrant_type,demtype) = dist_r(loc1+1,loc2+1,t1+1,immigrant_type,demtype) + surv(t1)*w2*w1*tr_z(j2,j4)*dist_hold;
                            end

                        end

                    end
                end
            end
        end
        
        if (amnesty~=0)&&(t1<Tr)
            amnesty_dist = squeeze(amnesty.*dist_1(:,:,:,t1+1,3,demtype)); 
            amnesty_dist = permute(amnesty_dist,[2,1,3]);
            amnesty_dist = squeeze(sum(amnesty_dist)); % amnesty_dist has dimensions nk,nb

            prod_legal = squeeze(dist_1(:,:,:,t1+1,2,demtype));   % has dimensions nk,nz,nb
            prod_legal = squeeze(sum(sum(permute(prod_legal,[1,3,2]))));
            prod_legal = prod_legal./sum(prod_legal);

            for i1 = 1:nk
                for i2 = 1:nb
                    for j2 = 1:nz
                        dist_1(i1,j2,i2,t1+1,2,demtype) = dist_1(i1,j2,i2,t1+1,2,demtype) + prod_legal(j2)*amnesty_dist(i1,i2);
                    end
                end
            end
            dist_1(:,:,:,t1+1,3,demtype) = (1-amnesty).*dist_1(:,:,:,t1+1,3,demtype);
        elseif (amnesty~=0)&&(t1==Tr)
            dist_r(:,:,t1+1,2,demtype) = dist_r(:,:,t1+1,2,demtype) + amnesty.*dist_r(:,:,t1+1,3,demtype);
            dist_r(:,:,t1+1,3,demtype) = (1-amnesty).*dist_r(:,:,t1+1,3,demtype);
        end

        dist_1(:,:,:,t1+1,3,demtype) = (1-deportation).*dist_1(:,:,:,t1+1,3,demtype);  
        dist_r(:,:,t1+1,3,demtype) = (1-deportation).*dist_r(:,:,t1+1,3,demtype);
        
    end


    % looping over dist_r
    if Tss-(tt1-1)>Tr
        for t1 = Tr+1:min(T-1,Tss-(tt1-1))
            %----------------------------------CHECK--------------------------
            %----------------------------------CHECK--------------------------
            %----------------------------------CHECK--------------------------
            age = t1;
            year = max(1,min(cohort +age - T,Tss)) +1;
            %----------------------------------CHECK--------------------------
            %----------------------------------CHECK--------------------------
            %----------------------------------CHECK--------------------------
%             im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(year); pop_trans(year)*imm_age(age)*illegal_rate(year)];   % using period 1 imm rate values for steady state
            im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(1); pop_trans(year)*imm_age(age)*illegal_rate(1)];   % using period 1 imm rate values for steady state
            for i1 = 1:nk
                for i2 = 1:nb
                    point_k = max(koptss(i1,i2,t1-Tr),kgrid(1));    % placing floor at lowest gridpoint.
                    loc1 = find(kgrid(1:nk-1)<=point_k,1,'last');
                    w1 = (point_k - kgrid(loc1))/(kgrid(loc1+1)-kgrid(loc1));  % amount allocated to higher gridpoint
                    w1 = min(w1,1);
                    for immigrant_type = 1:3
                        dist_hold = dist_r(i1,i2,t1,immigrant_type,demtype) + eq(i1,1)*eq(i2,1)*im_flow(immigrant_type);
                        
                        dist_r(loc1,i2,t1+1,immigrant_type,demtype) = dist_r(loc1,i2,t1+1,immigrant_type,demtype) + surv(t1)*(1-w1)*dist_hold;
                        dist_r(loc1+1,i2,t1+1,immigrant_type,demtype) = dist_r(loc1+1,i2,t1+1,immigrant_type,demtype) + surv(t1)*w1*dist_hold;
                    end
                end
            end
            if amnesty~=0 
                dist_r(:,:,t1+1,2,demtype) = dist_r(:,:,t1+1,2,demtype) + amnesty.*dist_r(:,:,t1+1,3,demtype);
                dist_r(:,:,t1+1,3,demtype) = (1-amnesty).*dist_r(:,:,t1+1,3,demtype);
            end
            
            dist_r(:,:,t1+1,3,demtype) = (1-deportation).*dist_r(:,:,t1+1,3,demtype);
        end
    end
    
    Kalive=zeros(1,min(Tss-(tt1-1),T));
    Kdead=zeros(1,min(Tss-(tt1-1),T));
    Lab = zeros(1,min(Tss-(tt1-1),T));
    ELab = zeros(1,min(Tss-(tt1-1),T));
    Dist= zeros(1,min(Tss-(tt1-1),T));
    Fedit = zeros(1,min(Tss-(tt1-1),T));
    SSrev = zeros(1,min(Tss-(tt1-1),T));
    SSexp = zeros(1,min(Tss-(tt1-1),T));
    Lfp = zeros(1,min(Tss-(tt1-1),T));
    SS_base = zeros(1,min(Tss-(tt1-1),T));



    % solving for aggregates by age
    for k1 = 1:min(Tss-(tt1-1),Tr)
        for i3 = 1:3
            for k2 = 1:nz
                for k4 = 1:nk
                    for k5 = 1:nb
                        Lab(k1) =Lab(k1)+ labopt(k4,k2,k5,k1)*dist_1(k4,k2,k5,k1,i3,demtype);
                        Kalive(k1) = Kalive(k1) + kopt(k4,k2,k5,k1)*dist_1(k4,k2,k5,k1,i3,demtype);
                        Kdead(k1) = Kdead(k1) + (1-surv(k1))*kopt(k4,k2,k5,k1)*dist_1(k4,k2,k5,k1,i3,demtype);
                        ELab(k1) =ELab(k1)+ z(k2,k1,demtype)*labopt(k4,k2,k5,k1)*dist_1(k4,k2,k5,k1,i3,demtype);
                        Dist(k1) = Dist(k1) + dist_1(k4,k2,k5,k1,i3,demtype);
                        Fedit(k1) = Fedit(k1)+ fitax(k4,k2,k5,k1)*dist_1(k4,k2,k5,k1,i3,demtype);
                        SSrev(k1) = SSrev(k1)+ fsstax(k4,k2,k5,k1)*dist_1(k4,k2,k5,k1,i3,demtype);
                        SSexp(k1) = SSexp(k1)+ max(ss_select(k4,k2,k5,k1)-1,0).*benopt(k4,k2,k5,k1)*dist_1(k4,k2,k5,k1,i3,demtype);
                        Lfp(k1) = Lfp(k1) + gt(labopt(k4,k2,k5,k1),0)*dist_1(k4,k2,k5,k1,i3,demtype);
                        SS_base(k1) = SS_base(k1)+ ss_base(k4,k2,k5,k1)*dist_1(k4,k2,k5,k1,i3,demtype);
                    end
                end
            end
        end
    end

    for k1 = Tr+1:min(T,Tss-(tt1-1))
        for i3 = 1:3
            for k4 = 1:nk
                for k5 = 1:nb
                    Kalive(k1) = Kalive(k1) + koptss(k4,k5,k1-Tr)*dist_r(k4,k5,k1,i3,demtype);
                    Kdead(k1) = Kdead(k1) + (1-surv(k1))*(koptss(k4,k5,k1-Tr)*dist_r(k4,k5,k1,i3,demtype));
                    Dist(k1) = Dist(k1) + dist_r(k4,k5,k1,i3,demtype);
                    Fedit(k1) = Fedit(k1)+ fitaxss(k4,k5,k1-Tr)*dist_r(k4,k5,k1,i3,demtype);
                    SSrev(k1) = SSrev(k1)+ fsstaxss(k4,k5,k1-Tr)*dist_r(k4,k5,k1,i3,demtype);
                    SSexp(k1) = SSexp(k1)+ benoptss(k4,k5,k1-Tr)*dist_r(k4,k5,k1,i3,demtype);
                end
            end
        end
    end
        
    filename = ['transvars_' num2str(demtype) '_' num2str(tt1) '_head_' num2str(polno) '.mat'];
    jobdir = ['job_' num2str(polno) '_' num2str(impolno)];  % gives name
    totfile = fullfile(jobdir,filename);
    save(totfile,'dist_1','dist_r','Kalive','Kdead','ELab','Lab','Dist','Fedit','SSrev','SSexp','Lfp','SS_base');
end


















