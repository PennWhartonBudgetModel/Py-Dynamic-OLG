function pop_trans = population_projector(impolno, jobdir)

load('params.mat');

load('Imm_Data.mat')

filename = ['imm_polparams_' num2str(impolno) '.mat'];
load(fullfile(jobdir, filename))



Tr = NRA_baseline;

load SSVALS pop_prev

demtype = 1;  % by the symmetry of the demographic type sizes, population will grow equally on each "island" (i.e., demtype island).  WLOG, we choose demtype=1.
filename = ['sspol' num2str(demtype) '.mat'];
totfile = fullfile(jobdir,filename);
load(totfile);

filename = ['distvars_' num2str(demtype) '.mat'];
totfile = fullfile(jobdir,filename);
load(totfile,'dist1ss','dist1ss_r');

dist1ss_previous = dist1ss;   % initiating working-age distribution array
dist1ss_r_previous = dist1ss_r;   % initiating retiree distribution array

pop = sum(sum(sum(sum(sum(dist1ss))))) + sum(sum(sum(sum(dist1ss_r))));    % initiating population

pop_trans = zeros(Tss,1);    % fill in with transition population sizes
pop_trans(1) = pop;

% dist_age_previous = ones(T,1);



for trans_year = 2:Tss
    
    fprintf('%u\n', trans_year)';
    
    year = trans_year - 1;
    
    dist1ss = zeros(nk,nz,nb,T,3);   % last dimension is [native, legal, illegal]
    dist1ss_r = zeros(nk,nb,T,3);
    
    % Initiating Natives
%             k0 = 0;    % initial capital  **********CHANGE IMMIGRANT CAPITAL IF THIS IS DIFFERENT FROM ZERO**********
%     im_flow = [pop_trans(year)*(pgr); pop_trans(year)*imm_age(1)*legal_rate(year); pop_trans(year)*imm_age(1)*illegal_rate(year)];   % using period 1 imm rate values for steady state
    im_flow = [pop_trans(year)*(pgr); pop_trans(year)*imm_age(1)*legal_rate(1); pop_trans(year)*imm_age(1)*illegal_rate(1)];   % using period 1 imm rate values for steady state

    loc1 = 1; loc2 = 1;    % initiating 
    for i1=1:nz
        dist1ss(loc1,i1,loc2,1,:) = squeeze(proddist_age(i1,1,:)).*(im_flow);
    end

    

    for t1 = 1:Tr
        age=t1;
%         im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(year); pop_trans(year)*imm_age(age)*illegal_rate(year)];
        im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(1); pop_trans(year)*imm_age(age)*illegal_rate(1)];
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
                                % lack of generality in im flow indicator below
                                dist_hold = squeeze(dist1ss_previous(i1,j2,i2,t1,immigrant_type)) + eq(i1,1)*eq(i2,1)*squeeze(proddist_age(j2,t1,immigrant_type)).*im_flow(immigrant_type);
                                
                                dist1ss(loc1,j4,loc2,t1+1,immigrant_type) = squeeze(dist1ss(loc1,j4,loc2,t1+1,immigrant_type)) + (surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                dist1ss(loc1+1,j4,loc2,t1+1,immigrant_type) = squeeze(dist1ss(loc1+1,j4,loc2,t1+1,immigrant_type)) + (surv(t1)*(1-w2)*(w1)*tr_z(j2,j4)).*dist_hold;
                                dist1ss(loc1,j4,loc2+1,t1+1,immigrant_type) = squeeze(dist1ss(loc1,j4,loc2+1,t1+1,immigrant_type)) + (surv(t1)*(w2)*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                dist1ss(loc1+1,j4,loc2+1,t1+1,immigrant_type) = squeeze(dist1ss(loc1+1,j4,loc2+1,t1+1,immigrant_type)) + (surv(t1)*(w2)*(w1)*tr_z(j2,j4)).*dist_hold;
                            end
                        elseif t1==Tr
                            for immigrant_type = 1:3
                                dist_hold = squeeze(dist1ss_previous(i1,j2,i2,t1,immigrant_type)) + eq(i1,1)*eq(i2,1)*squeeze(proddist_age(j2,t1,immigrant_type)).*im_flow(immigrant_type);
                                
                                dist1ss_r(loc1,loc2,t1+1,immigrant_type) = squeeze(dist1ss_r(loc1,loc2,t1+1,immigrant_type)) + (surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                dist1ss_r(loc1+1,loc2,t1+1,immigrant_type) = squeeze(dist1ss_r(loc1+1,loc2,t1+1,immigrant_type)) + (surv(t1)*(1-w2)*(w1)*tr_z(j2,j4)).*dist_hold;
                                dist1ss_r(loc1,loc2+1,t1+1,immigrant_type) = squeeze(dist1ss_r(loc1,loc2+1,t1+1,immigrant_type)) + (surv(t1)*(w2)*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                dist1ss_r(loc1+1,loc2+1,t1+1,immigrant_type) = squeeze(dist1ss_r(loc1+1,loc2+1,t1+1,immigrant_type)) + (surv(t1)*(w2)*(w1)*tr_z(j2,j4)).*dist_hold;
                            end
                        end
                    end
                end
            end
        end
    end

    for t1 = Tr+1:T-1
        age=t1;
%         im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(year); pop_trans(year)*imm_age(age)*illegal_rate(year)];
        im_flow = [0; pop_trans(year)*imm_age(age)*legal_rate(1); pop_trans(year)*imm_age(age)*illegal_rate(1)];
        for i1 = 1:nk
            for i2 = 1:nb
                point_k = max(koptss(i1,i2,t1-Tr),kgrid(1));    % placing floor at lowest gridpoint.
                loc1 = find(kgrid(1:nk-1)<=point_k,1,'last');
                w1 = (point_k - kgrid(loc1))/(kgrid(loc1+1)-kgrid(loc1));  % amount allocated to higher gridpoint
                w1 = min(w1,1);
                for immigrant_type = 1:3
                    dist_hold = squeeze(dist1ss_r_previous(i1,i2,t1,immigrant_type)) + eq(i1,1)*eq(i2,1)*im_flow(immigrant_type);
                    
                    dist1ss_r(loc1,i2,t1+1,immigrant_type) = squeeze(dist1ss_r(loc1,i2,t1+1,immigrant_type)) + (surv(t1)*(1-w1)).*dist_hold;
                    dist1ss_r(loc1+1,i2,t1+1,immigrant_type) = squeeze(dist1ss_r(loc1+1,i2,t1+1,immigrant_type)) + (surv(t1).*w1)*dist_hold;
                end
            end
        end
    end
    
    dist1ss(:,:,:,:,3) = (1-deportation).*dist1ss(:,:,:,:,3);  
    dist1ss_r(:,:,:,3) = (1-deportation).*dist1ss_r(:,:,:,3);


    pop = sum(sum(sum(sum(sum(dist1ss))))) + sum(sum(sum(sum(dist1ss_r))));    % summing population
%             dist1ss = dist1ss.*(T/dist_tot);    % normalizing the distribution so that it has size T
%             dist1ss_r = dist1ss_r.*(T/dist_tot);
    pop_trans(trans_year) = pop;

%     dist_age = [squeeze(sum(sum(sum(sum(permute(dist1ss,[1,2,3,5,4]))))))+squeeze(sum(sum(sum(permute(dist1ss_r,[1,2,4,3])))))];
%     death_rate = sum(dist_age'.*(1-surv(1:T)))/pop
%     disteps = max(abs((1/dist_age(1)).*dist_age(2:end) - (1/dist_age_previous(1)).*dist_age_previous(2:end)))
%     dist_age_previous = dist_age;
%             disteps = max(max(max(max(max(abs(dist1ss_previous - dist1ss)))))) + max(max(max(max(abs(dist1ss_r_previous - dist1ss_r)))))
    dist1ss_previous = dist1ss;   % updating working-age distribution array
    dist1ss_r_previous = dist1ss_r;   % updating retiree distribution array

end

% pop_trans = [pop_prev; pop_trans];