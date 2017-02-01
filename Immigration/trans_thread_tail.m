function trans_thread_tail(tt1,polno,impolno)
% tt1=47;
% polno=1;
% impolno=1;


load params.mat
% load Pop_Evol.mat
load Surv_Probs.mat
load Imm_Data.mat
% Alive at start of transition


jobdir = 'Testing';
load(fullfile(jobdir, sprintf('imm_polparams_%u.mat', impolno)))


% amnesty=.2;
hello=2;
MU2 = zeros(ndem,T);
% MU_ILL = zeros(ndem,T);
% MU_LEG = zeros(ndem,T);
for demtype = 1:ndem
    MU2(demtype,:) = (demdist_2015(demtype)).*((1/sum(Mu2)).*Mu2);    % going to repeat last value twice
%     MU_ILL(demtype,:) = ((1/sum(Mu2)).*Mu2);
%     MU_LEG(demtype,:) = ((1/sum(Mu2)).*Mu2);
end
% plot(mu2)
surv = 1-surv_proj(1,:);
MU3 = zeros(ndem,T);
for t1 = 1:T
    for i1 = 1:ndem
        MU3(i1,t1) = (1-surv(t1))*MU2(i1,t1);
    end
end


% immigration
% if length(im_rate)<Tss
%     leg_rate = [leg_rate',leg_rate(end).*ones(1,T)];
% end

% leg_rate2 = leg_rate; %[im_rate(tt1:Tss), im_rate(end).*ones(1,T-(Tss-(tt1-1)))];
% if length(im_rate)<Tss
%     illegal_rate = [illegal_rate',illegal_rate(end).*ones(1,T)];
% end

% illeg_rate2 = illegal_rate; %[im_rate(tt1:Tss), im_rate(end).*ones(1,T-(Tss-(tt1-1)))];
% imm_age2 = [imm_age(tt1:Tss), imm_age(end).*ones(1,T-(Tss-(tt1-1)))];
% pgrowth2 = pgrowth; %[pgrowth(tt1:Tss), pgrowth(end).*ones(1,T-(Tss-(tt1-1)))];

% dyi = dist_year_ill(1,:);
% dyl = dist_year_leg(1,:);


% load DIST_DATA1_debug
% jobdir = ['job_' num2str(1) '_' num2str(1)];
filename = ['eqmdist'];
totfile = fullfile(jobdir,filename);
load(totfile);


% tt1
dist_1 = zeros(nk,nz,nb,T-(tt1-1)+1,3,ndem);
dist_r = zeros(nk,nb,T-(tt1-1)+1,3,ndem);

cohort = T-(tt1-1);
load SSVALS pop_prev
pop_trans = [pop_prev; pop_trans];
% legal_rate = [legal_rate(1); legal_rate];
% illegal_rate = [illegal_rate(1); illegal_rate];

for demtype = 1:ndem
%     demtype=2;
    mu2 = MU2(demtype,:);
    mu3 = MU3(demtype,:);
    
%     mi2 = MU_ILL(demtype,:);
%     ml2 = MU_LEG(demtype,:);

    filename = ['tail' num2str(tt1) '_' num2str(demtype) '_' num2str(polno)  '.mat'];
    totfile = fullfile('Freeze', 'Cohorts', filename);
    load(totfile);
    
    
    % while beqeps>beqtol
    
    dist_1(:,:,:,1,1,demtype) = dist1(:,:,:,tt1,1,demtype);  % Natives
    dist_1(:,:,:,1,2,demtype) = dist1(:,:,:,tt1,2,demtype);  % Legals
    dist_1(:,:,:,1,3,demtype) = dist1(:,:,:,tt1,3,demtype);  % Illegals
    
    dist_r(:,:,1,1,demtype) = distr(:,:,tt1,1,demtype);   % Natives
    dist_r(:,:,1,2,demtype) = distr(:,:,tt1,2,demtype);   % Legals
    dist_r(:,:,1,3,demtype) = distr(:,:,tt1,3,demtype);   % Illegals
    

    if tt1<=Tr
        for t1 = tt1:Tr
            %----------------------------------CHECK--------------------------
            %----------------------------------CHECK--------------------------
            %----------------------------------CHECK--------------------------
            age = t1;
            year = max(1,min(cohort +age - T,Tss)) +1;
            %----------------------------------CHECK--------------------------
            %----------------------------------CHECK--------------------------
            %----------------------------------CHECK--------------------------
%             year = cohort +age - T;
%             im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(year); pop_trans(year)*imm_age(age)*illegal_rate(year)];
            im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(1); pop_trans(year)*imm_age(age)*illegal_rate(1)];
            node=0;
            for j2 = 1:nz
                for i1 = 1:nk
                    node=node+1;
                    for i2 = 1:nb
                        
                        point_k = max(kopt(i1,j2,i2,t1 - (tt1-1)),kgrid(1));    % placing floor at lowest gridpoint.
                        loc1 = find(kgrid(1:nk-1)<=point_k,1,'last');
                        w1 = (point_k - kgrid(loc1))/(kgrid(loc1+1)-kgrid(loc1));  % amount allocated to higher gridpoint
                        w1 = min(w1,1);

                        point_b = max(bopt(i1,j2,i2,t1 - (tt1-1)),bgrid(1));    % placing floor at lowest gridpoint.
                        loc2 = find(bgrid(1:nb-1)<=point_b,1,'last');    % lower gridpoint
                        w2 = (point_b - bgrid(loc2))/(bgrid(loc2+1)-bgrid(loc2));  % amount allocated to higher gridpoint
                        w2 = min(w2,1);
                        
                        for j4=1:nz
                            if t1<Tr
                                for immigrant_type = 1:3
                                    dist_hold = dist_1(i1,j2,i2,t1 - (tt1-1),immigrant_type,demtype) + eq(i1,1)*eq(i2,1)*proddist_age(j2,t1,immigrant_type)*im_flow(immigrant_type);
                                    
                                    dist_1(loc1,j4,loc2,t1 - (tt1-1) + 1,immigrant_type,demtype) = dist_1(loc1,j4,loc2,t1 - (tt1-1) + 1,immigrant_type,demtype) + surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)*dist_hold;
                                    dist_1(loc1+1,j4,loc2,t1 - (tt1-1) + 1,immigrant_type,demtype) = dist_1(loc1+1,j4,loc2,t1 - (tt1-1) + 1,immigrant_type,demtype) + surv(t1)*(1-w2)*(w1)*tr_z(j2,j4)*dist_hold;
                                    dist_1(loc1,j4,loc2+1,t1 - (tt1-1) + 1,immigrant_type,demtype) = dist_1(loc1,j4,loc2+1,t1 - (tt1-1) + 1,immigrant_type,demtype) + surv(t1)*w2*(1-w1)*tr_z(j2,j4)*dist_hold;
                                    dist_1(loc1+1,j4,loc2+1,t1 - (tt1-1) + 1,immigrant_type,demtype) = dist_1(loc1+1,j4,loc2+1,t1 - (tt1-1) + 1,immigrant_type,demtype) + surv(t1)*w2*(w1)*tr_z(j2,j4)*dist_hold;
                                end

                            elseif t1==Tr
                                for immigrant_type = 1:3
                                    dist_hold = dist_1(i1,j2,i2,t1 - (tt1-1),immigrant_type,demtype) + eq(i1,1)*eq(i2,1)*proddist_age(j2,t1,immigrant_type)*im_flow(immigrant_type);
                                    
                                    dist_r(loc1,loc2,t1 - (tt1-1) + 1,immigrant_type,demtype) = dist_r(loc1,loc2,t1 - (tt1-1) + 1,immigrant_type,demtype) + surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)*dist_hold;
                                    dist_r(loc1+1,loc2,t1 - (tt1-1) + 1,immigrant_type,demtype) = dist_r(loc1+1,loc2,t1 - (tt1-1) + 1,immigrant_type,demtype) + surv(t1)*(1-w2)*w1*tr_z(j2,j4)*dist_hold;
                                    dist_r(loc1,loc2+1,t1 - (tt1-1) + 1,immigrant_type,demtype) = dist_r(loc1,loc2+1,t1 - (tt1-1) + 1,immigrant_type,demtype) + surv(t1)*w2*(1-w1)*tr_z(j2,j4)*dist_hold;
                                    dist_r(loc1+1,loc2+1,t1 - (tt1-1) + 1,immigrant_type,demtype) = dist_r(loc1+1,loc2+1,t1 - (tt1-1) + 1,immigrant_type,demtype) + surv(t1)*w2*w1*tr_z(j2,j4)*dist_hold;
                                end

                            end
                            
                        end
                    end
                end
            end
            if (amnesty~=0)&&(t1<Tr)
                amnesty_dist = squeeze(amnesty.*dist_1(:,:,:,t1 - (tt1-1) + 1,3,demtype)); 
                amnesty_dist = permute(amnesty_dist,[2,1,3]);
                amnesty_dist = squeeze(sum(amnesty_dist)); % amnesty_dist has dimensions nk,nb
                
                prod_legal = squeeze(dist_1(:,:,:,t1 - (tt1-1) +1,2,demtype));   % has dimensions nk,nz,nb
                prod_legal = squeeze(sum(sum(permute(prod_legal,[1,3,2]))));
                prod_legal = prod_legal./sum(prod_legal);
                
                for i1 = 1:nk
                    for i2 = 1:nb
                        for j2 = 1:nz
                            dist_1(i1,j2,i2,t1 - (tt1-1) +1,2,demtype) = dist_1(i1,j2,i2,t1 - (tt1-1) +1,2,demtype) + prod_legal(j2)*amnesty_dist(i1,i2);
                        end
                    end
                end
                dist_1(:,:,:,t1 - (tt1-1) + 1,3,demtype) = (1-amnesty).*dist_1(:,:,:,t1 - (tt1-1) + 1,3,demtype);
            elseif (amnesty~=0)&&(t1==Tr)
                dist_r(:,:,t1 - (tt1-1) + 1,2,demtype) = dist_r(:,:,t1 - (tt1-1) + 1,2,demtype) + amnesty.*dist_r(:,:,t1 - (tt1-1) + 1,3,demtype);
                dist_r(:,:,t1 - (tt1-1) + 1,3,demtype) = (1-amnesty).*dist_r(:,:,t1 - (tt1-1) + 1,3,demtype);
            end
            
            dist_1(:,:,:,t1 - (tt1-1) + 1,3,demtype) = (1-deportation).*dist_1(:,:,:,t1 - (tt1-1) + 1,3,demtype);  
            dist_r(:,:,t1 - (tt1-1) + 1,3,demtype) = (1-deportation).*dist_r(:,:,t1 - (tt1-1) + 1,3,demtype);
        end
    end


    % looping over dist_r
    for t1 = max(Tr+1,tt1):T-1
        %----------------------------------CHECK--------------------------
        %----------------------------------CHECK--------------------------
        %----------------------------------CHECK--------------------------
        age = t1;
        year = max(1,min(cohort +age - T,Tss))+1;
        %----------------------------------CHECK--------------------------
        %----------------------------------CHECK--------------------------
        %----------------------------------CHECK--------------------------
%         year = cohort +age - T;
%         im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(year); pop_trans(year)*imm_age(age)*illegal_rate(year)];
        im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(1); pop_trans(year)*imm_age(age)*illegal_rate(1)];
        for i1 = 1:nk
            for i2 = 1:nb
                point_k = max(koptss(i1,i2,t1-max(Tr,tt1-1)),kgrid(1));    % placing floor at lowest gridpoint.
                loc1 = find(kgrid(1:nk-1)<=point_k,1,'last');
                w1 = (point_k - kgrid(loc1))/(kgrid(loc1+1)-kgrid(loc1));  % amount allocated to higher gridpoint
                w1 = min(w1,1);
                for immigrant_type = 1:3
                    dist_hold = dist_r(i1,i2,t1 - (tt1-1),immigrant_type,demtype) + eq(i1,1)*eq(i2,1)*im_flow(immigrant_type);
                    
                    dist_r(loc1,i2,t1 - (tt1-1)+1,immigrant_type,demtype) = dist_r(loc1,i2,t1 - (tt1-1)+1,immigrant_type,demtype) + surv(t1)*(1-w1)*dist_hold;
                    dist_r(loc1+1,i2,t1 - (tt1-1)+1,immigrant_type,demtype) = dist_r(loc1+1,i2,t1 - (tt1-1)+1,immigrant_type,demtype) + surv(t1)*w1*dist_hold; 
                end
            end
        end
        if amnesty~=0 
            dist_r(:,:,t1 - (tt1-1) + 1,3,demtype) = (1-amnesty).*dist_r(:,:,t1 - (tt1-1) + 1,3,demtype);
            dist_r(:,:,t1 - (tt1-1) + 1,2,demtype) = dist_r(:,:,t1 - (tt1-1) + 1,2,demtype) + amnesty.*dist_r(:,:,t1 - (tt1-1) + 1,3,demtype);
        end
        dist_r(:,:,t1 - (tt1-1) + 1,3,demtype) = (1-deportation).*dist_r(:,:,t1 - (tt1-1) + 1,3,demtype);
    end

%         Kalive=zeros(1,min(Tss-(tt1-1),T));
%         Kdead=zeros(1,min(Tss-(tt1-1),T));
%         Lab = zeros(1,min(Tss-(tt1-1),T));
%         ELab = zeros(1,min(Tss-(tt1-1),T));
%         Dist= zeros(1,min(Tss-(tt1-1),T));
    Kalive=zeros(1,T);
    Kdead=zeros(1,T);
    Lab = zeros(1,T);
    ELab = zeros(1,T);
    Dist= zeros(1,T);
    Fedit = zeros(1,T);
    SSrev = zeros(1,T);
    SSexp = zeros(1,T);
    Lfp = zeros(1,T);
    SS_base = zeros(1,T);
%     SSrev = zeros(1,T);


%         for k1 = tt1:min(Tss,T)
%             node1 = 1;
%             for k2 = 1:nz
%                 for k4 = 1:nk
%                     for k5 = 1:nb
%                         Kalive(k1-(tt1-1)) = Kalive(k1-(tt1-1)) + mu2(k1)*kgrid(k4)*(dist_1(node1,k5,k1 - (tt1-1),demtype)+dist_r(node1,k5,k1 - (tt1-1),demtype));
%                         Kdead(k1-(tt1-1)) = Kdead(k1-(tt1-1)) + mu3(k1)*kgrid(k4)*(dist_1(node1,k5,k1 - (tt1-1),demtype)+dist_r(node1,k5,k1 - (tt1-1),demtype));
%                     end
%                     node1 = node1+1;
%                 end
%             end
%         end


%         KPR(tt1+Tss,1:min(Tss,T)-(tt1-1)) = KPR(tt1+Tss,1:min(Tss,T)-(tt1-1))+Kalive(1:min(Tss,T)-(tt1-1))+Kdead(1:min(Tss,T)-(tt1-1));
%         BEQ(tt1+Tss,1:min(Tss,T)-(tt1-1)) = BEQ(tt1+Tss,1:min(Tss,T)-(tt1-1))+Kdead(1:min(Tss,T)-(tt1-1));



    % solving for aggregates by age: region 1
    if tt1<=Tr
        for k1 = tt1:Tr
            for i3 = 1:3
                for k2 = 1:nz
                    for i1 = 1:nk
                        node=node+1;
                        for i2 = 1:nb
                            Lab(k1-(tt1-1)) =Lab(k1-(tt1-1))+ labopt(i1,k2,i2,k1-(tt1-1))*dist_1(i1,k2,i2,k1-(tt1-1),i3,demtype);
                            Kalive(k1-(tt1-1)) = Kalive(k1-(tt1-1)) + kopt(i1,k2,i2,k1-(tt1-1))*dist_1(i1,k2,i2,k1-(tt1-1),i3,demtype);
                            Kdead(k1-(tt1-1)) = Kdead(k1-(tt1-1)) + (1-surv(k1))*kopt(i1,k2,i2,k1-(tt1-1))*dist_1(i1,k2,i2,k1-(tt1-1),i3,demtype);
                            ELab(k1-(tt1-1)) =ELab(k1-(tt1-1))+ z(k2,k1,demtype)*labopt(i1,k2,i2,k1-(tt1-1))*dist_1(i1,k2,i2,k1-(tt1-1),i3,demtype);
                            Dist(k1-(tt1-1)) = Dist(k1-(tt1-1)) + dist_1(i1,k2,i2,k1-(tt1-1),i3,demtype);
                            Fedit(k1-(tt1-1)) = Fedit(k1-(tt1-1))+ fitax(i1,k2,i2,k1-(tt1-1))*dist_1(i1,k2,i2,k1-(tt1-1),i3,demtype);
                            SSrev(k1-(tt1-1)) = SSrev(k1-(tt1-1))+ fsstax(i1,k2,i2,k1-(tt1-1))*dist_1(i1,k2,i2,k1-(tt1-1),i3,demtype);
                            SSexp(k1-(tt1-1)) = SSexp(k1-(tt1-1))+ benopt(i1,k2,i2,k1-(tt1-1))*dist_1(i1,k2,i2,k1-(tt1-1),i3,demtype);
                            Lfp(k1-(tt1-1)) =Lfp(k1-(tt1-1))+ gt(labopt(i1,k2,i2,k1-(tt1-1)),0)*dist_1(i1,k2,i2,k1-(tt1-1),i3,demtype);
                            SS_base(k1-(tt1-1)) = SS_base(k1-(tt1-1))+ ss_base(i1,k2,i2,k1-(tt1-1))*dist_1(i1,k2,i2,k1-(tt1-1),i3,demtype);
                        end
                    end
                end
            end
        end
    end

    % solving for aggregates by age: region 3
%     T-(tt1-1)
    for k1 = max(Tr+1,tt1):T%min(T,Tss)
        for i3 = 1:3
            for i1 = 1:nk
                for i2 = 1:nb
                    Kalive(k1-(tt1-1)) = Kalive(k1-(tt1-1)) + koptss(i1,i2,k1-max(Tr,tt1-1))*dist_r(i1,i2,k1-(tt1-1),i3,demtype);
                    Kdead(k1-(tt1-1)) = Kdead(k1-(tt1-1)) + (1-surv(k1))*(koptss(i1,i2,k1-max(Tr,tt1-1))*dist_r(i1,i2,k1-(tt1-1),i3,demtype));
                    Dist(k1-(tt1-1)) = Dist(k1-(tt1-1)) + dist_r(i1,i2,k1-(tt1-1),i3,demtype);
                    Fedit(k1-(tt1-1)) = Fedit(k1-(tt1-1))+ fitaxss(i1,i2,k1-max(Tr,tt1-1))*dist_r(i1,i2,k1-(tt1-1),i3,demtype);
                    SSrev(k1-(tt1-1)) = SSrev(k1-(tt1-1))+ fsstaxss(i1,i2,k1-max(Tr,tt1-1))*dist_r(i1,i2,k1-(tt1-1),i3,demtype);
                    SSexp(k1-(tt1-1)) = SSexp(k1-(tt1-1))+ benoptss(i1,i2,k1-max(Tr,tt1-1))*dist_r(i1,i2,k1-(tt1-1),i3,demtype);
                end
            end
        end
    end
    
    filename = ['transvars_' num2str(demtype) '_' num2str(tt1) '_tail_' num2str(polno) '.mat'];
    totfile = fullfile(jobdir,filename);
    save(totfile,'dist_1','dist_r','Kalive','Kdead','ELab','Lab','Dist','Fedit','SSrev','SSexp','Lfp','SS_base');

end






