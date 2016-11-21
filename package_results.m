%%
% Package all results into csv files for staging.
% 
%%


function [] = package_results()

% Identify csv save directory
csv_dir = dirFinder.csv;

% Clear or create csv directory
if exist(csv_dir, 'dir')
    rmdir(csv_dir, 's')
end
mkdir(csv_dir)


% Load target elasticity sets
s = load(fullfile(dirFinder.param, 'ss_inverses.mat'));
elasticity_sets = s.targets(:,[2,3]);
deep_param_sets = s.inverses;
clear('s')


% Specify identifier strings and numerical codes for aggregates
agg_codes = { 102, 'feditlab'         ;
              103, 'ssrev'            ;
              110, 'fcaprev'          ;
              111, 'domestic_fcaptax' ;
              112, 'foreign_fcaptax'  ;
              401, 'domestic_cap'     ;
              402, 'foreign_cap'      ;
              403, 'domestic_debt'    ;
              404, 'foreign_debt'     ;
              199, 'Y'                ;
              202, 'ssexp'            ;
              299, 'Y'                ;
                1, 'cap'              ;
                6, 'Y'                ;
               24, 'elab'             ;
               25, 'lfpr'             ;
               28, 'labinc'           ;
               29, 'kinc'             };
n_aggs = size(agg_codes, 1);

% Specify projection years for csv files
years_csv = (2015 : 2089)';
n_years = length(years_csv);

% Specify number of years to shift results
% (Currently, shifting up from 2015 to 2017)
upshift = 2;


% Load table of run IDs
tablefilename = 'scenario_table.csv';
tablesource = fullfile(dirFinder.root, 'charts', 'scenario_table', tablefilename);
copyfile(tablesource, fullfile(csv_dir, tablefilename))

fid = fopen(tablesource);
idtable = textscan(fid, ['%u64 ', repmat('%*q ', 1, 7),'%q %*[^\n]'], 'Delimiter', ',', 'HeaderLines', 1);
fclose(fid);

n_ids = length(idtable{1});

% Define mapping from table plan names to dynamic model plan names
planmap = struct( 'CurrentLaw'  , 'base'    , ...
                  'Clinton'     , 'clinton' , ...
                  'Trump'       , 'trump'   , ...
                  'Ryan'        , 'ryan'    );


for i = 1:n_ids
    
    % Get run ID
    id = idtable{1}(i);
    
    % Get run specifications
    specs = xml2struct(idtable{2}{i});
    specs = specs.Dynamics;
    
    % Skip run if not dynamic model tax policy
    if ~isfield(specs, 'TaxPlan'), continue, end
    
    % Get economy openness
    openness = str2double(specs.OpenEconomy.Text);
    
    % Skip run if not fully open or fully closed
    if (openness ~= 1 && openness ~= 0), continue, end
    
    % Get elasticities
    labor_elas   = str2double(specs.LaborElasticity  .Text);
    savings_elas = str2double(specs.SavingsElasticity.Text);
    
    % Find elasticity set index, which is equivalent to the deep parameter set index
    inddeep = all(bsxfun(@eq, [labor_elas, savings_elas], elasticity_sets), 2);
    
    % Get deep parameters based on set index
    deep_params = deep_param_sets(inddeep,:);
    
    % Get plan
    plan = planmap.(specs.TaxPlan.Text);
    
    % Get government expenditure reduction
    % (Note negation to align with dynamic model convention)
    gcut = -str2double(specs.ExpenditureShift.Text);
    
    % Get dynamic baseline flag
    % (Note activation for closed economy runs only)
    use_dynamic_baseline = str2double(specs.UseDynamicBaseline.Text) & (openness == 0);
    
    
    % Identify working directories
    switch openness
        case 1, save_dir = dirFinder.open  (deep_params, plan, gcut);
        case 0, save_dir = dirFinder.closed(deep_params, plan, gcut);
    end
    
    
    % Identify iterations log file
    logfile = fullfile(save_dir, 'iterations.txt');
    
    % Check for run completion
    if exist(logfile, 'file')
        
        % Get last line of iterations log file
        fid = fopen(logfile);
        while ~feof(fid)
            lastline = fgetl(fid);
        end
        fclose(fid);
        
        % Check convergence error against threshold
        lastiter = sscanf(lastline, '  %2d  --  %f', [1,2]);
        if (lastiter(2) > 0.1), continue, end
        
    else
        continue
    end
    
    
    % Load aggregates    
    s_dynamic = load(fullfile(save_dir, 'aggregates.mat'       ));
    s_static  = load(fullfile(save_dir, 'aggregates_static.mat'));
    
    if (use_dynamic_baseline)
        s_base_open   = load(fullfile(dirFinder.open  (deep_params, 'base', +0.00), 'aggregates.mat'));
        s_base_closed = load(fullfile(dirFinder.closed(deep_params, 'base', +0.00), 'aggregates.mat'));
    end
    
    
    % Find number of projection years
    Tss = length(s_dynamic.kpr_total);
    
    % Find number of entries to be trimmed or padded
    trim_or_pad = upshift + Tss - n_years;
    n_trim =  max(trim_or_pad, 0);
    n_pad  = -min(trim_or_pad, 0);
    
    % Extract aggregates and generate csvs
    for j = 1:n_aggs
        
        aggnum = agg_codes{j,1};
        aggstr = agg_codes{j,2};
        
        agg_dynamic = s_dynamic.([aggstr,'_total' ]);
        agg_static  = s_static .([aggstr,'_static']);
        
        % Adjust static aggregate if using dynamic baseline
        if (use_dynamic_baseline)
            agg_static  = agg_static .* s_base_open.([aggstr,'_total']) ./ s_base_closed.([aggstr,'_total']);
            agg_static(isnan(agg_static)) = 0;
        end
        
        agg_series  = [ ones(upshift, 2); ...
                        [agg_dynamic(1:end-n_trim)', agg_static(1:end-n_trim)']; ...
                        ones(n_pad,   2) ];
        
        csvfile = fullfile(csv_dir, sprintf('%u-%u.csv', id, aggnum));
        fid = fopen(csvfile, 'w');
        fprintf(fid, 'Year,DynamicAggregate,StaticAggregate\n');
        fclose(fid);
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


fprintf('\nResults successfully packaged into csv files:\n\t%s\n', csv_dir)


end




%% xml2struct
%
% Pulled from Matlab File Exchange on 2016-11-02:
% 
%   https://www.mathworks.com/matlabcentral/fileexchange/28518
%   https://www.mathworks.com/matlabcentral/fileexchange/58700
% 
% (Note use of Java for string inputs)
%

function [outStruct] = xml2struct(input)
% input can be a Java XML object, an XML file, or a string in
% XML format.
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
% 
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Originally written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increase by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
% Modified by X. Mo, University of Wisconsin, 12-5-2012
% Modified by Chao-Yuan Yeh, August 2016

errorMsg = ['%s is not in a supported format.\n\nInput has to be',...
        ' a java xml object, an xml file, or a string in xml format.'];

% check if input is a java xml object
if isa(input, 'org.apache.xerces.dom.DeferredDocumentImpl') ||...
        isa(input, 'org.apache.xerces.dom.DeferredElementImpl')
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

% parse xDoc into a MATLAB structure
outStruct = parseChildNodes(xDoc);
    
end

function [children, ptext, textflag] = parseChildNodes(theNode)
% Recurse over node children.
children = struct;
ptext = struct; 
textflag = 'Text';

if hasChildNodes(theNode)
    childNodes = getChildNodes(theNode);
    numChildNodes = getLength(childNodes);

    for count = 1:numChildNodes

        theChild = item(childNodes,count-1);
        [text, name, attr, childs, textflag] = getNodeData(theChild);
        
        if ~strcmp(name,'#text') && ~strcmp(name,'#comment') && ...
                ~strcmp(name,'#cdata_dash_section')
            % XML allows the same elements to be defined multiple times,
            % put each in a different cell
            if (isfield(children,name))
                if (~iscell(children.(name)))
                    % put existsing element into cell format
                    children.(name) = {children.(name)};
                end
                index = length(children.(name))+1;
                % add new element
                children.(name){index} = childs;
                
                textfields = fieldnames(text);
                if ~isempty(textfields)
                    for ii = 1:length(textfields)
                        children.(name){index}.(textfields{ii}) = ...
                            text.(textfields{ii});
                    end
                end
                if(~isempty(attr)) 
                    children.(name){index}.('Attributes') = attr; 
                end
            else
                % add previously unknown (new) element to the structure
                
                children.(name) = childs;
                
                % add text data ( ptext returned by child node )
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

            % this is the text in an element (i.e., the parentNode) 
            if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                    ptext.(ptextflag) = text.(textflag);
                else
                    % This is what happens when document is like this:
                    % <element>Text <!--Comment--> More text</element>
                    %
                    % text will be appended to existing ptext
                    ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                end
            end
        end

    end
end
end

function [text,name,attr,childs,textflag] = getNodeData(theNode)
% Create structure of node info.

%make sure name is allowed as structure name
name = toCharArray(getNodeName(theNode))';
name = strrep(name, '-', '_dash_');
name = strrep(name, ':', '_colon_');
name = strrep(name, '.', '_dot_');
name = strrep(name, '_', 'u_');

attr = parseAttributes(theNode);
if (isempty(fieldnames(attr))) 
    attr = []; 
end

%parse child nodes
[childs, text, textflag] = parseChildNodes(theNode);

% Get data from any childless nodes. This version is faster than below.
if isempty(fieldnames(childs)) && isempty(fieldnames(text))
    text.(textflag) = toCharArray(getTextContent(theNode))';
end

% This alterative to the above 'if' block will also work but very slowly.
% if any(strcmp(methods(theNode),'getData'))
%   text.(textflag) = char(getData(theNode));
% end
    
end

function attributes = parseAttributes(theNode)
% Create attributes structure.
attributes = struct;
if hasAttributes(theNode)
   theAttributes = getAttributes(theNode);
   numAttributes = getLength(theAttributes);

   for count = 1:numAttributes
        % Suggestion of Adrian Wanner
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



