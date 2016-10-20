clear
close('all')
clc

addpath('..')

% Get deep parameters corresponding to specified linear index


plans = {'base','clinton','trump','ryan'};


% Specify expenditure reduction
gcuts ={ -0.05, +0.00, +0.05, +0.10};


% Specify open/closed
econ_cases = {'open','closed'};

var_cases = { 'cap_total'         ,  'Total Capital'                  ;
              'elab_total'        ,  'Efficient Labor'                ;
              'Y_total'           ,  'Output'                         ;
              'kpr_total'         ,  'Aggregate Saving'               ;
              'lab_total'         ,  'Labor Hours'                    ;
              'beq_total'         ,  'Aggregate Bequests'             ;
              'lfpr_total'        ,  'Labor Force Participation Rate' ;
              'debt_total'        ,  'Total Debt'                     ;
              'fedit_total'       ,  'Personal Income Tax Revenue'    ;
              'ssrev_total'       ,  'Payroll Tax Revenue'            ;
              'fcaptax_total'     ,  'Capital Tax Revenue'            ;
              'ssexp_total'       ,  'Social Security Expenditures'   ;
              'Gtilde'            ,  'G-tilde'                        ;
              'Ttilde'            ,  'T-tilde'                        ;
              'fedincome_total'   ,  'Personal Income Tax Base'       ;
              'feditlab_total'    ,  'PIT Tax Revenue from Labor'     ;
              'labinc_total'      ,  'Total Labor Income'             ;
              'fcaprev_total'     ,  'Tax Revenue from Capital'       ;
              'kinc_total'        ,  'Total Capital Income'           };
             
             


% open only values
var_cases_open = { 'domestic_cap_total'     ,  'Domestic Capital'               ; 
                   'domestic_debt_total'    ,  'Domestic Debt'                  ;
                   'domestic_fcaptax_total' ,  'Domestic Capital Tax Revenue'   ;
                   'foreign_cap_total'      ,  'Foreign Capital'                ;
                   'foreign_debt_total'     ,  'Foreign Debt'                   ;
                   'foreign_fcaptax_total'  ,  'Foreign Capital Tax Revenue'    };
              
% closed only values
var_cases_closed = { 'rhos'      ,   'Capital-to-Labor Ratio'       ;
                     'wages'     ,   'Wages'                        ;
                     'cap_shares',   'Capital Share of Portfolio'   ;
                     'rate_caps' ,   'Capital Return (including q)' };

nvars = length(var_cases(:,1)) ;             
              
s = load(fullfile('..','Parameters','param_global.mat'));
Tss = s.Tss;
time = 2016+1:1:2016+Tss;
clear('s')


for inddeep = 6%1:16
    
    deep_parameters = inddeep_to_params(inddeep);
    for econ_case = econ_cases

        % create vectors for and case variables and names 
        if strcmp(econ_case,'open')
            var_cases_total = [var_cases(:,1)', var_cases_open(:,1)'];
            names           = [var_cases(:,2)', var_cases_open(:,2)'];
        elseif strcmp(econ_case,'closed')
            var_cases_total = [var_cases(:,1)', var_cases_closed(:,1)'];
            names           = [var_cases(:,2)', var_cases_closed(:,2)'];
        end

        for gcut = gcuts

            % Change figure file name here for separate files along "gcut", delete if exists
            if strcmp(econ_case,'open')
                figure_pdf = [sprintf('params=%d',inddeep) '_open.pdf']; 
            elseif strcmp(econ_case,'closed')
                figure_pdf = [sprintf('params=%d',inddeep) '_' sprintf('closed_gcut=%+0.2f', gcut{1}) '.pdf'];
            end

            if exist(figure_pdf) %#ok<EXIST>
                delete(figure_pdf)
            end

            % Compiling data corresponding to variable, econ_case, plan and creating graphic
            for var_case = var_cases_total
                var_plot = [];

                for plan = plans 
                    % Find save directory for open economy transition path and load aggregates
                    if     strcmp(econ_case{1}, 'open'  )
                        load_dir     = dirFinder.open  (deep_parameters(1), deep_parameters(2), deep_parameters(3), plan{1});
                        figure_title = 'Open Economy';
                    elseif strcmp(econ_case{1}, 'closed')
                        load_dir     = dirFinder.closed(deep_parameters(1), deep_parameters(2), deep_parameters(3), plan{1}, gcut{1});
                        figure_title = ['Closed Economy, Government Expenditure Cut = ' num2str(gcut{1}*100) '%'];
                    end

                    load_aggregates   = load(fullfile(load_dir, 'aggregates.mat'));
                    if strcmp(econ_case,'closed') && ~isempty(find(strcmp(var_cases_closed(:,1),var_case{1}))) %#ok<EFIND>
                        load_aggregates   = load(fullfile(load_dir, 'solution.mat'));
                    end
                    var_plot = [var_plot; load_aggregates.(var_case{1})]; %#ok<AGROW>
                end   


                fig1 = figure(1);
                plot(time,var_plot,'LineWidth',4)
                title(figure_title)
                ylabel(names{find(strcmp(var_cases_total,var_case{1}))}); %#ok<FNDSB>
                legend(plans,'Location','Northwest','Orientation','horizontal')
                fig1.PaperUnits = 'normalized';
                fig1.PaperPosition = [.1,.1,.9,.9];
                fig1.PaperOrientation = 'landscape';
                print(fig1,'figure1','-dpdf')
                append_pdfs(figure_pdf,'figure1.pdf')

            end

            % Exit loop after first iteration if open economy
            if strcmp(econ_case{1},'open')
                break
            end
        end
    end
end

if exist('figure1.pdf') %#ok<EXIST>
    delete('figure1.pdf')
end



