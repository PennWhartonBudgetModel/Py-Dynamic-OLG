%%
% Dynamic model production run manager.
% 
%%


classdef modelProducer

properties (Constant)
    
    % Define run directory and run file path
    run_dir  = fullfile(dirFinder.source(), 'Runs');
    run_file = @(irun) fullfile(modelProducer.run_dir, sprintf('run%04d.mat', irun));
    
end

methods (Static)
    
    % Define full set of production runs
    function [] = define_runs()
        
        % Clear or create run directory
        if exist(modelProducer.run_dir, 'dir'), rmdir(modelProducer.run_dir, 's'), end, mkdir(modelProducer.run_dir)
        
        % Construct elasticity inverter
        [~, invert] = modelCalibrator.invert();
        
        % Define vectors of elasticities, to be inverted into baseline parameters
        labelasv = num2cell(0.00 : 0.20 : 1.00);
        savelasv = num2cell(0.00 : 0.20 : 1.00);
        
        % Define vectors of counterfactual parameters
        taxplanv = {'trumpA', 'trumpB'};
        gcutv    = num2cell(+0.20 : -0.01 : 0.00);
        
        % Generate all parameter combinations
        [labelasinds, savelasinds, taxplaninds, gcutinds] = ndgrid(1:length(labelasv), 1:length(savelasv), 1:length(taxplanv), 1:length(gcutv));
        
        % Construct baseline and counterfactual definitions for runs and save
        for irun = 1:numel(labelasinds)
            
            basedef    = invert(struct('labelas', labelasv{labelasinds(irun)}, 'savelas', savelasv{savelasinds(irun)})); %#ok<NASGU>
            counterdef = struct('taxplan', taxplanv{taxplaninds(irun)}, 'gcut', gcutv{gcutinds(irun)});                  %#ok<NASGU>
            
            save(modelProducer.run_file(irun), 'basedef', 'counterdef');
            
        end
        
    end
    
    
    % Perform production run
    function [] = run(irun)
        
        % Load run baseline and counterfactual definitions
        s = load(modelProducer.run_file(irun));
        
        % Execute dynamic model solver
        % (Closed economy will generate corresponding steady state and open economy dependencies)
        save_dir = dynamicSolver.closed(s.basedef, s.counterdef);
        
        % Extract and save solver termination condition
        iterations = csvread(fullfile(save_dir, 'iterations.csv'));
        
        termination.iter = iterations(end, 1);
        termination.eps  = iterations(end, 2);
        
        termination = termination; %#ok<ASGSL,NASGU>
        save(modelProducer.run_file(irun), '-append', 'termination');
        
    end
    
    
    % Check production run termination conditions
    function [] = check_terminations()
        
        % Identify production run definition files
        run_files = dir(fullfile(modelProducer.run_dir, 'run*.mat'));
        nrun = length(run_files);
        
        % Initialize cell array of termination conditions
        terminations = cell(0,5);
        
        for irun = 1:nrun
            
            % Load run definition and termination condition
            s = load(modelProducer.run_file(irun));
            
            % Generate baseline and counterfactual definition tags
            [basedef_tag, counterdef_tag] = dynamicSolver.generate_tags(s.basedef, s.counterdef);
            
            % Store termination condition, setting default values if missing
            if ~isfield(s, 'termination'), s.termination = struct('iter', Inf, 'eps', Inf); end
            terminations = [terminations; {irun, basedef_tag, counterdef_tag, s.termination.iter, s.termination.eps}]; %#ok<AGROW>
            
        end
        
        % Sort termination conditions by increasing iterations and error terms
        [~, sortinds] = sortrows(cell2mat(terminations(:,4:5)));
        
        % Save termination conditions to csv file
        fid = fopen(fullfile(dirFinder.saveroot(), 'terminations.csv'), 'w');
        fprintf(fid, 'Run,Baseline Definition,Counterfactual Definition,Termination Iteration,Termination Error Term\n');
        for irun = 1:nrun, fprintf(fid, '%d,%s,%s,%d,%0.4f\n', terminations{sortinds(irun),:}); end
        fclose(fid);
        
    end
    
    
    % Package production run results into csv files for front end deployment
    function [] = package_results()
        
        % Identify csv directory
        csv_dir = dirFinder.csv();
    
        % Clear or create csv directory
        if exist(csv_dir, 'dir'), rmdir(csv_dir, 's'), end, mkdir(csv_dir)
    
    
        % Specify identifier strings and numerical codes for aggregates
        codes = { 102, 'labpits'        ;
                  103, 'ssts'           ;
                  110, 'caprevs'        ;
                  111, 'cits_domestic'  ;
                  112, 'cits_foreign'   ;
                  401, 'caps_domestic'  ;
                  402, 'caps_foreign'   ;
                  403, 'debts_domestic' ;
                  404, 'debts_foreign'  ;
                  199, 'outs'           ;
                  202, 'bens'           ;
                  299, 'outs'           ;
                    1, 'caps'           ;
                    6, 'outs'           ;
                   24, 'labeffs'        ;
                   25, 'lfprs'          ;
                   28, 'labincs'        ;
                   29, 'capincs'        };
        
        % Specify projection years for csv files
        years_csv = (2016 : 2089)';
        nyear = length(years_csv);
        
        % Specify number of years to shift results
        % (Currently, shifting up from 2016 to 2018)
        nshift = 2;
        
                
        % Load table of run IDs
        tablefile = 'scenario_table.csv';
        tablesource = fullfile(dirFinder.root, 'charts', 'scenario_table', tablefile);
        copyfile(tablesource, fullfile(csv_dir, tablefile))
        
        fid = fopen(tablesource);
        idtable = textscan(fid, ['%u64 ', repmat('%*q ', 1, 7), '%q %*[^\n]'], 'Delimiter', ',', 'HeaderLines', 1);
        fclose(fid);
        
        nid = length(idtable{1});
        
        % Define mapping from table tax plan names to dynamic model names
        taxplanmap = struct( 'CurrentLaw'  , 'base'    , ...
                             'TrumpA'      , 'trumpA'  , ...
                             'TrumpB'      , 'trumpB'  );
        
        % Construct elasticity inverter
        [~, invert] = modelCalibrator.invert();
        
        % Initialize missing aggregates flag
        missing = false;
        
        
        for i = 1:nid
            
            fprintf('Processing run ID %6d of %6d\n', i, nid)
            
            % Get run ID
            id = idtable{1}(i);
            
            % Get run definition
            def = xml2struct(idtable{2}{i});
            def = def.Dynamics;
            
            % Skip run if tax plan absent or not currently modeled
            if ~isfield(def, 'TaxPlan') || ~isfield(taxplanmap, def.TaxPlan.Text), continue, end
            
            % Get economy openness and skip run if not fully open or fully closed
            switch str2double(def.OpenEconomy.Text)
                case 1, economy = 'open'  ;
                case 0, economy = 'closed';
                otherwise, continue
            end
            
            % Get elasticities
            labelas = str2double(def.LaborElasticity  .Text);
            savelas = str2double(def.SavingsElasticity.Text);
            
            % Invert elasticities to get baseline definition
            basedef = invert(struct('labelas', labelas, 'savelas', savelas));
            
            % Get tax plan
            taxplan = taxplanmap.(def.TaxPlan.Text);
            isbase = strcmp(taxplan, 'base');
            
            % Get government expenditure reduction
            % (Note negation to align with dynamic model convention and enforcement of positive zero)
            gcut = -str2double(def.ExpenditureShift.Text);
            if (gcut == 0), gcut = +0.00; end
            
            % Construct counterfactual definition
            if isbase
                counterdef = struct();
            else
                counterdef = struct('taxplan', taxplan, 'gcut', gcut);
            end
            
            % Get dynamic baseline flag
            % (Applicable to closed economy runs only)
            usedynamicbaseline = str2double(def.UseDynamicBaseline.Text) & strcmp(economy, 'closed');
            
            
            % Identify working directories
            save_dir = dirFinder.save(economy, basedef, counterdef);
            
            % Load aggregates
            try
                Dynamic = load(fullfile(save_dir, 'dynamics.mat'));
                if ~isbase
                    Static = load(fullfile(save_dir, 'statics.mat'));
                end
            catch
                missing = true;
                continue
            end
            
            if usedynamicbaseline
                Dynamic_open_base   = load(fullfile(dirFinder.save('open'  , basedef), 'dynamics.mat'));
                Dynamic_closed_base = load(fullfile(dirFinder.save('closed', basedef), 'dynamics.mat'));
            end
            
            
            % Find number of model years
            T_model = length(Dynamic.(codes{1,2}));
            
            % Find number of entries to be trimmed or padded
            nextra = nshift + T_model - nyear;
            ntrim  =  max(nextra, 0);
            npad   = -min(nextra, 0);
            
            
            for j = 1:size(codes,1)
                
                % Get aggregate name
                name = codes{j,2};
                
                % Extract dynamic and static aggregate series
                a.dynamic = Dynamic.(name);
                if isbase
                    a.static = a.dynamic;
                else
                    a.static = Static.(name);
                end
                
                % Adjust static aggregate series if using dynamic baseline
                if usedynamicbaseline
                    a.static = a.static .* Dynamic_open_base.(name) ./ Dynamic_closed_base.(name);
                    a.static(isnan(a.static)) = 0;
                end
                
                % Consolidate, shift, trim, and pad aggregate series
                agg_series  = [ ones(nshift, 2); 
                                [a.dynamic(1:end-ntrim)', a.static(1:end-ntrim)'];
                                ones(npad,   2) ];
                
                % Save aggregate series to csv file
                csvfile = fullfile(csv_dir, sprintf('%u-%u.csv', id, codes{j,1}));
                fid = fopen(csvfile, 'w'); fprintf(fid, 'Year,DynamicAggregate,StaticAggregate\n'); fclose(fid);
                dlmwrite(csvfile, [years_csv, agg_series], '-append')
                
            end
            
        end
        
        
        % Get time stamp from csv directory name and convert format
        [~, timestamp] = fileparts(csv_dir);
        timestamp = datestr(datenum(timestamp, 'yyyy-mm-dd-HH-MM'), 'yyyy-mm-dd HH:MM');
        
        % Compose tag file
        fid = fopen(fullfile(csv_dir, 'dynamicModelTag.txt'), 'w');
        fprintf(fid, 'TimeStamp=%s\nShortDescription=Generated by dynamic model version %s.', timestamp, dirFinder.get_commit_id);
        fclose(fid);
        
        
        fprintf('\nResults packaged into csv files:\n\t%s\n', csv_dir)
        if missing, warning('Some aggregates not found'), end
        
    end % package_results
    
    
    %% 
    % Small run of immigration counter-factuals
    function [] = immigration_run()
        
         % Construct elasticity inverter
         [~, invert] = modelCalibrator.invert();
         % Invert elasticities to get baseline definition, using "default"
         basedef = invert(struct('labelas', 0.5, 'savelas', 0.5));
         
         % Calculate prem_legal from requested policies
         prod_skilled   = 1.4;
         current_legal  = 0.45;
         prod_unskilled = 37/55;  
         % from prod_skilled*current_legal + prod_unskilled*(1-current_legal) = prod_immigrants 
         % rem: in baseline, prod_immigrants assumed = 1
         
         for immScale = [0.6, 0.5]
            for portion_skilled = [0.55, 0.75]
                immPremium = prod_skilled*portion_skilled + prod_unskilled*(1-portion_skilled);
                counterdef = struct(    'legal_scale'   , immScale    ...
                                    ,   'prem_legal'    , immPremium  ...
                                    );
                fprintf( '--------------------------------\n' );
                fprintf( 'RUNNING immScale=%f, portion_skilled=%f \n ', immScale, portion_skilled );
                dynamicSolver.closed( basedef, counterdef, '' );               
            end
         end
         
    end % immigration run
    
end % methods

end % modelProducer




%% xml2struct
%
% Pulled from Matlab File Exchange on 2016-11-02:
% 
%   https://www.mathworks.com/matlabcentral/fileexchange/28518
%   https://www.mathworks.com/matlabcentral/fileexchange/58700
% 
% (Note use of Java for string inputs)
%

function [output] = xml2struct(input)
% Input can be a Java XML object, an XML file, or a string in XML format.
% 
% XML:
% 
%   <XMLname attrib1="Some value">
%     <Element>Some text</Element>
%     <DifferentElement attrib2="2">Some more text</Element>
%     <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
%   </XMLname>
% 
% Matlab structure:
% 
%   s.XMLname.Attributes.attrib1 = "Some value";
%   s.XMLname.Element.Text = "Some text";
%   s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
%   s.XMLname.DifferentElement{1}.Text = "Some more text";
%   s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
%   s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
%   s.XMLname.DifferentElement{2}.Text = "Even more text";
% 
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
% 
% Originally written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increase by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
% Modified by X. Mo, University of Wisconsin, 12-5-2012
% Modified by Chao-Yuan Yeh, August 2016

errorMsg = '%s is not in a supported format.\n\nInput has to be a java xml object, an xml file, or a string in xml format.';
if isa(input, 'org.apache.xerces.dom.DeferredDocumentImpl') || isa(input, 'org.apache.xerces.dom.DeferredElementImpl')
    xDoc = input;
else
    try 
        if exist(input, 'file') == 2
            xDoc = xmlread(input);
        else
            try
                xDoc = xmlFromString(input);
            catch
                error(errorMsg, inputname(1));
            end
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
            error(errorMsg, inputname(1));
        else
            rethrow(ME)
        end
    end
end
output = parseChildNodes(xDoc);
end

function [children, ptext, textflag] = parseChildNodes(theNode)
children = struct;
ptext = struct; 
textflag = 'Text';
if hasChildNodes(theNode)
    childNodes = getChildNodes(theNode);
    numChildNodes = getLength(childNodes);
    for count = 1:numChildNodes
        theChild = item(childNodes,count-1);
        [text, name, attr, childs, textflag] = getNodeData(theChild);
        if ~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section')
            if (isfield(children,name))
                if (~iscell(children.(name)))
                    children.(name) = {children.(name)};
                end
                index = length(children.(name))+1;
                children.(name){index} = childs;
                textfields = fieldnames(text);
                if ~isempty(textfields)
                    for ii = 1:length(textfields)
                        children.(name){index}.(textfields{ii}) = text.(textfields{ii});
                    end
                end
                if(~isempty(attr)) 
                    children.(name){index}.('Attributes') = attr; 
                end
            else
                children.(name) = childs;
                textfields = fieldnames(text);
                if ~isempty(textfields)
                    for ii = 1:length(textfields)
                        children.(name).(textfields{ii}) = text.(textfields{ii});
                    end
                end
                if(~isempty(attr)) 
                    children.(name).('Attributes') = attr; 
                end
            end
        else
            ptextflag = 'Text';
            if (strcmp(name, '#cdata_dash_section'))
                ptextflag = 'CDATA';
            elseif (strcmp(name, '#comment'))
                ptextflag = 'Comment';
            end
            if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                    ptext.(ptextflag) = text.(textflag);
                else
                    ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                end
            end
        end
    end
end
end

function [text,name,attr,childs,textflag] = getNodeData(theNode)
name = toCharArray(getNodeName(theNode))';
name = strrep(name, '-', '_dash_');
name = strrep(name, ':', '_colon_');
name = strrep(name, '.', '_dot_');
name = strrep(name, '_', 'u_');
attr = parseAttributes(theNode);
if (isempty(fieldnames(attr))) 
    attr = []; 
end
[childs, text, textflag] = parseChildNodes(theNode);
if isempty(fieldnames(childs)) && isempty(fieldnames(text))
    text.(textflag) = toCharArray(getTextContent(theNode))';
end
end

function attributes = parseAttributes(theNode)
attributes = struct;
if hasAttributes(theNode)
   theAttributes = getAttributes(theNode);
   numAttributes = getLength(theAttributes);
   for count = 1:numAttributes
        str = toCharArray(toString(item(theAttributes,count-1)))';
        k = strfind(str,'='); 
        attr_name = str(1:(k(1)-1));
        attr_name = strrep(attr_name, '-', '_dash_');
        attr_name = strrep(attr_name, ':', '_colon_');
        attr_name = strrep(attr_name, '.', '_dot_');
        attributes.(attr_name) = str((k(1)+2):(end-1));
   end
end
end

function xmlroot = xmlFromString(iString)
import org.xml.sax.InputSource
import java.io.*
iSource = InputSource();
iSource.setCharacterStream(StringReader(iString));
xmlroot = xmlread(iSource);
end