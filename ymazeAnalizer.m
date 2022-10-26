function [Summary,Stats,Stats2] = ymazeAnalizer(varargin)
% Y-Maze data analyzer: Analizes data from the Y-watermaze for reversal
% learning in mice
%
% Max van der Doe, 2022
% ErasmusUniversity, Doctor Molewaterplein 40, 3015 GD Rotterdam
% m.vanderdoe.1@erasmusmc.nl
%
%Function = [Summary, Stats2] = ymazeAnalizer(varargin)
%_______________________________________________________________________________
% Input: all input arguments must be paired, order does not matter 
%        e.g. ymazeAnalizer('File','[File name..]','Condition1','Genotype','SupressPlot','True')
%
% File:             provide name of file to be analysed. file must be an excel file
%                   Providing file name as input argument nesesetates the
%                   fle be added to matlab path.
%
% Condition1/2/3:   Provide the condition by which the analisis wil be
%                   done. the maximum anount of conditions is 3.
%                   add 1,2,3 after the dondition to analize more conditions.
%                   e.g. 'Condition1','Condition2','Condition3'
%
% Condition inputs: 'Genotype', 'Housing','Sex'
%
% Supressplot:      Supress the generation of plots. use this when only the
%                   statistical data is of relevents, or the programe is
%                   used within another programe. input must be true or
%                   false
%________________________________________________________________________________
% Output:
%
% Summary:          All meaned values based on group and condition(s)
%                   provided. Data includes group count, mean correct choices for each sesion,
%                   mean value where the mice reach a succes rate of 80%
%                   (mean_Acq80 % mean_Rev80) and the mean area under the
%                   curve (mean_AUCAcq & mean_AUCRev)
% 
% Stats2            All values for each individual mouse. Data includes
%                   Gene, Gender and the provided condition(s), mean values
%                   for each session, session number where the succes rate
%                   was more than 80% (Acq80 % Rev80) and the Area under
%                   the curve (AUCAcq &AUCRev)
%
% Saveing Data...   Output data is also saved to the provided file. Summary and stats2 are 
%                   saved to the 'Summary' and 'Statisticks' sheets respectivly. 
%                   when these sheets do not exist they are created
%________________________________________________________________________________
% dependencies       shadedErrorBar.m 
%                    Rob Campbell (2022). raacampbell/shadedErrorBar 
%                    (https://github.com/raacampbell/shadedErrorBar), GitHub. 
%                    Retrieved March 4, 2022.
%                    https://nl.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
%
%                    violin.m
%                    Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel
%                    density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
%                    hhoffmann@uni-bonn.de
%                    https://nl.mathworks.com/matlabcentral/fileexchange/45134-violin-plot


%% check function input arguments
Condition = {};                                                             %PreDefine Condition to be empty
ConditionOptions = {'Genotype','Housing','Sex'};                            %predefine Condition options
SupressPlot = [];                                                           %perdefine SupressPlot to be empty

%resolve File name
if isempty(find(strcmpi(varargin,'File'),1))==0                             % If file name is given, than:
    File = varargin{find(strcmpi(varargin,'File'))+1};                      % File = the input argument after 'File'
    try
        isfile(File);                                                       % Check if the provided file is a valid file
    catch
        warning('Input must be a valid file')                               % If not, the user points to one
        File= uigetfile;
    end
    Path = fileparts(which(File));                                          % get path of file
    File = [Path '\' File];                                                 % Combine path and file, this is used to save the data

else
    [File,Path] = uigetfile('*');                                                % Get file from user
    File = [Path File];                                                     % Combine path and file, this is used to save the data
end

%Resolve condition
if isempty(find(strcmpi(varargin,'Condition1'),1))==0                       % If Condition1 exists than:
    Condition{1} = varargin{find(strcmpi(varargin,'Condition1'))+1};        % Condition{1} = the input directly after Condition1
    if isempty(find(strcmpi(varargin,'Condition2'),1))==0                   % If Condition2 exists than:
        Condition{2} = varargin{find(strcmpi(varargin,'Condition2'))+1};    % Condition{2} = the input directly after Condition2
        if isempty(find(strcmpi(varargin,'Condition3'),1))==0               % If Condition3 exists than:
            Condition{3} = varargin{find(strcmpi(varargin,'Condition3'))+1};% Condition{3} = the input directly after Condition3
        end
    end
end

%if condition is not resolved, ask for condition as a list
if isempty(Condition)                                                       % Condition is predefined as empty, therefor if no input arguments are provided this wil be true
    disp('please select a condition...')
    [indx, tf]= listdlg('PromptString',{'Select a condition.','multiple conditions can be selected',''},'ListString',ConditionOptions,'SelectionMode','multiple'); % Ask for user imput
    if tf == 0, error('No Condition was selected, please try again..'), end %if the user pressed cancel, give error
    Condition(1) = ConditionOptions(indx(1));                               % resolve condition same as above
    if length(indx)>=2
        Condition(2) = ConditionOptions(indx(2));
        if length(indx)>=3
            Condition(3) = ConditionOptions(indx(3));
        end
    end
end

%resolve plot supression
if isempty(find(strcmpi(varargin,'SupressPlot'),1))==0                      % If SupressPlot is given as an imput argument:
    if strcmpi(varargin{find(strcmpi(varargin,'SupressPlot'))+1},'true')    % 'true' = true
        SupressPlot = true; 
    elseif strcmpi(varargin{find(strcmpi(varargin,'SupressPlot'))+1},'false')% 'false' = false
        SupressPlot = false;
    else                                                                    % If input argument is not recocnised ask user for input
        SupressPlot = questdlg('would you like to supress the plots?','Plot Supression','Yes','No','No');
        switch SupressPlot                                                  % switch based on response
            case 'Yes'                                                      % 'yes' = true
                SupressPlot = true;
            case 'No'
                SupressPlot = false;                                        % 'no' = false
        end
    end
end
if isempty(SupressPlot)                                                     % SupressPlot is predefined as empty, therefor if no input arguments are provided this wil be true
    SupressPlot = questdlg('would you like to supress the plots?','Plot Supression','Yes','No','No');% ask for user input
    switch SupressPlot                                                      %switch based on response
        case 'Yes'                                                          % 'yes' = true
            SupressPlot = true;
        case 'No'                                                           % 'no' = false
            SupressPlot = false;
    end
end
%% import the datasheet
disp('Loading data...')
data = readtable(File,"UseExcel",true);                                     % inport data as a table
data = data(4:(sum(find(any(strcmp(data.MouseNumber,''),2),4))-7),(1:195)); % trim table to include only the relevant data

%fill al the blank cells left after importing merged cell
for i = 1:height(data)                                                      % loop over the first row
    if strcmp(data.Group(i),'')                                             % see if  the row is empty
        data.Group(i) = data.Group(i-1);                                    % if so than replace '' by the value of the previous row
    end
end

%% Determine the condition of analysis
Stats = array2table(zeros(height(data.MouseNumber),length(Condition)+25));  % Predefine Stats to be a table, width is determint by the number of conditions given
switch length(Condition)                                                    % Swich based on the amount of conditions
    case 1                                                                  % Variable names are the same as trings in the Condition variable
        Stats.Properties.VariableNames = {'MouseNumber','Gene',Condition{1},'Acq1.1','Acq1.2','Acq1.3','Acq1.4','Acq2.1','Acq2.2','Acq2.3','Acq2.4','Test','Rev1.1','Rev1.2','Rev1.3','Rev1.4','Rev1.5','Rev2.1','Rev2.2','Rev2.3','Rev2.4','Rev2.5','Acq80','Rev80','AUCAcq','AUCRev'};
    case 2
        Stats.Properties.VariableNames = {'MouseNumber','Gene',Condition{1},Condition{2},'Acq1.1','Acq1.2','Acq1.3','Acq1.4','Acq2.1','Acq2.2','Acq2.3','Acq2.4','Test','Rev1.1','Rev1.2','Rev1.3','Rev1.4','Rev1.5','Rev2.1','Rev2.2','Rev2.3','Rev2.4','Rev2.5','Acq80','Rev80','AUCAcq','AUCRev'};
    case 3
        Stats.Properties.VariableNames = {'MouseNumber','Gene',Condition{1},Condition{2},Condition{3},'Acq1.1','Acq1.2','Acq1.3','Acq1.4','Acq2.1','Acq2.2','Acq2.3','Acq2.4','Test','Rev1.1','Rev1.2','Rev1.3','Rev1.4','Rev1.5','Rev2.1','Rev2.2','Rev2.3','Rev2.4','Rev2.5','Acq80','Rev80','AUCAcq','AUCRev'};
end

% fill columns with the relevant information
Stats.MouseNumber = data.MouseNumber;
Stats.Gene = data.Group;                                                    % Import the Group from data, this is used to group the data
Stats.(Condition{1}) = table2cell(data(:,find(strcmpi(Condition{1},data.Properties.VariableNames))));
if length(Condition)>=2
    Stats.(Condition{2}) = table2cell(data(1:end,find(strcmpi(Condition{2},data.Properties.VariableNames))));
end
if length(Condition)>=3
    Stats.(Condition{3}) = table2cell(data(1:end,find(strcmpi(Condition{3},data.Properties.VariableNames))));
end

%% calculate mean values
disp('Calculating...')
start = find(strcmpi(Stats.Properties.VariableNames,'Acq1.1'));             % find the column where the data starts,'Acq1.1'. this is used to index

for y = 1:height(data.MouseNumber)                                          % Loop over every row
    for x = start:start+18                                                  % Loop over evely collumn.
        Stats(y,x) = {sum(data{y,(10*(x-start))+6:2:(10*(x-start))+15})/5}; % Calculate the mean for each sesion and place in the stats variable
    end

    if isempty(find(table2array(Stats(y,start:start+12))>=.8,1))
        Stats.Acq80(y) = 9;
    else
        Stats.Acq80(y) = find(table2array(Stats(y,start:start+12))>=.8,1);
    end
    if isempty(find(table2array(Stats(y,start+13:start+22))>=.8,1))
        Stats.Rev80(y) = 19;
    else
        Stats.Rev80(y) = find(table2array(Stats(y,start+13:start+22))>=.8,1);
    end
    
    if isnan(Stats{y,start+8})
        yAcq = [table2array(Stats(y,start:start+4)) table2array(Stats(y,start+8))];
    else
        yAcq = table2array(Stats(y,start:start+8));
    end
    yRev = table2array(Stats(y,start+9:start+18));

    Stats.AUCAcq(y) = trapz(yAcq);
    Stats.AUCRev(y) = trapz(yRev);

end

Stats2 = Stats(Stats.Test >= .8,:);                                          % exclude mice with a test score of lower than .8

switch length(Condition)                                                    % switch based on number of conditions
    case 1                                                                  % Calculate Mean by gene and condition
        Summary = groupsummary(Stats2,{'Gene',Condition{1}},'Mean',vartype("numeric")); 
    case 2
        Summary = groupsummary(Stats2,{'Gene',Condition{1},Condition{2}},'Mean',vartype("numeric"));
    case 3
        Summary = groupsummary(Stats2,{'Gene',Condition{1},Condition{2},Condition{3}},'Mean',vartype("numeric"));
end

%% test for NaN and if found replace with Mean_AUC from the Summary table

Rows = isnan(Stats2.AUCAcq);
Stats2.AUCAcq(Rows) = trapz(table2array(Stats2(Rows,start:start+3)),2);
Stats2.AUCAcq(isnan(Stats2.AUCAcq)) = 0;

Stats2.AUCRev(isnan(Stats2.AUCRev)) = 0;

%% Line Plot and violin plot statistics

switch SupressPlot
    case false
        Groups = unique(Summary.Gene);
        disp('Generating plots...')
        Row=1; %Row indicates the row that needs to be plotted from summary
        l=1;
        r=1;
        for i = 1:2:(height(unique(Summary.Gene))*2)
            figure(i)
            fprintf('Figure %d.\n',i)
            
            for j = 1:2:sum(strcmp(Summary.Gene,Groups(r)))*2
                if length(Condition)==1,Rows = (strcmp(Summary.Gene(Row),Stats2.Gene)&strcmp(Summary.(Condition{1})(Row),Stats2.(Condition{1}))); end
                if length(Condition)==2,Rows = (strcmp(Summary.Gene(Row),Stats2.Gene)&strcmp(Summary.(Condition{1})(Row),Stats2.(Condition{1}))&strcmp(Summary.(Condition{2})(Row),Stats2.(Condition{2}))); end
                if length(Condition)==3,Rows = (strcmp(Summary.Gene(Row),Stats2.Gene)&strcmp(Summary.(Condition{1})(Row),Stats2.(Condition{1}))&strcmp(Summary.(Condition{2})(Row),Stats2.(Condition{2}))&strcmp(Summary.(Condition{3})(Row),Stats2.(Condition{3}))); end

                StdAcq = std(table2array(Stats2(Rows,start:start+8)),0,1,'omitnan');
                StdRev = std(table2array(Stats2(Rows,start+9:start+18)),0,1,'omitnan');

                subplot(sum(strcmp(Summary.Gene,Groups(r))),2,j)
                fprintf('Subplot %d.\n',j)
                plot(1:9,table2array(Summary(Row,start:start+8)))
                title('Aquisition '+string(Summary.Gene(Row))+' '+strjoin(table2cell(Summary(Row,2:length(Condition)+1)))+' '+'n='+string(Summary.GroupCount(Row)))
                shadedErrorBar(1:9,table2array(Summary(Row,start:start+8)),StdAcq,'lineProps','b')
                xlim([1 9]),ylim([0 1])
                xline(Summary.mean_Acq80(Row),'-',{'More than','80% learned',Summary.mean_Acq80(Row)})

                subplot(sum(strcmp(Summary.Gene,Groups(r))),2,j+1)
                fprintf('Subplot %d.\n',j+1)
                plot(1:10,table2array(Summary(Row,start+9:start+18)))
                title('Reversal '+string(Summary.Gene(Row))+' '+strjoin(table2cell(Summary(Row,2:length(Condition)+1)))+' '+'n='+string(Summary.GroupCount(Row)))
                shadedErrorBar(1:10,table2array(Summary(Row,start+9:start+18)),StdRev,'lineProps','b')
                xlim([1 10]),ylim([0 1])
                xline(Summary.mean_Rev80(Row),'-',{'More than','80% learned',Summary.mean_Rev80(Row)})

                Row = Row+1;
            end
            figure(i+1)
            fprintf('Figure %d.\n',i+1)
            for j = 1:2:sum(strcmp(Summary.Gene,Groups(r)))*2

                if length(Condition)==1,Rows = (strcmp(Summary.Gene(l),Stats2.Gene)&strcmp(Summary.(Condition{1})(l),Stats2.(Condition{1}))); end
                if length(Condition)==2,Rows = (strcmp(Summary.Gene(l),Stats2.Gene)&strcmp(Summary.(Condition{1})(l),Stats2.(Condition{1}))&strcmp(Summary.(Condition{2})(l),Stats2.(Condition{2}))); end
                if length(Condition)==3,Rows = (strcmp(Summary.Gene(l),Stats2.Gene)&strcmp(Summary.(Condition{1})(l),Stats2.(Condition{1}))&strcmp(Summary.(Condition{2})(l),Stats2.(Condition{2}))&strcmp(Summary.(Condition{3})(l),Stats2.(Condition{3}))); end

                subplot(sum(strcmp(Summary.Gene,Groups(r))),2,j)
                fprintf('Subplot %d.\n',j)
                violin(Stats2.AUCAcq(Rows,:));
                title('Aquisition '+string(Summary.Gene(l))+' '+strjoin(table2cell(Summary(l,2:length(Condition)+1)))+' '+'n='+string(Summary.GroupCount(l)))
                ylim([1 10])
                scatter((rand(1,height(Stats2.AUCAcq(Rows,:)))*(1.25-.75)+.75),Stats2.AUCAcq(Rows,:),'red','filled','^')
                legend({'AUC','Mean','Median','Mice'})

                subplot(sum(strcmp(Summary.Gene,Groups(r))),2,j+1)
                fprintf('Subplot %d.\n',j+1)
                violin(Stats2.AUCRev(Rows,:));
                title('Reversal '+string(Summary.Gene(l))+' '+strjoin(table2cell(Summary(l,2:length(Condition)+1)))+' '+'n='+string(Summary.GroupCount(l)))
                ylim([1 10])
                scatter((rand(1,height(Stats2.AUCRev(Rows,:)))*(1.25-.75)+.75),Stats2.AUCRev(Rows,:),'red','filled','^')
                legend({'AUC','Mean','Median','Mice'})

                l=l+1;
            end
            r=r+1;
        end
    case true
        disp('Plots are supressed')
end

%% generate output file
fid = -1;
while fid==-1
    fid = fopen(File,'a');
    if fid==-1
        questdlg('Unable to save data. Check if file is open and press ok to continue','Saving File','ok');
    else
        fclose(fid);
        disp('Saveing data...')
        writetable(Summary,File,'UseExcel',true,'Sheet','Grouped summary','WriteVariableNames',true,'writemode','overwritesheet')
        writetable(Stats,File,'UseExcel',true,'Sheet','Statisticks every mouse','writemode','overwritesheet','WriteVariableNames',true)
        writetable(Stats2,File,'UseExcel',true,'Sheet','Statisticks more than 80%','writemode','overwritesheet','WriteVariableNames',true)
    end
end

msgbox('Analysis succesfull')
end