function tableout = importfile(workbookFile,sheetName,startRow,endRow)
%IMPORTFILE Import data from a spreadsheet
%   DATA = IMPORTFILE(FILE) reads data from the first worksheet in the
%   Microsoft Excel spreadsheet file named FILE and returns the data as a
%   table.
%
%   DATA = IMPORTFILE(FILE,SHEET) reads from the specified worksheet.
%
%   DATA = IMPORTFILE(FILE,SHEET,STARTROW,ENDROW) reads from the specified
%   worksheet for the specified row interval(s). Specify STARTROW and
%   ENDROW as a pair of scalars or vectors of matching size for
%   dis-contiguous row intervals. To read to the end of the file specify an
%   ENDROW of inf.%
% Example:
%   UOK262UMRC6HK2BioRx2peakscontroldata09152017 = importfile('UOK262_UMRC6_HK2_BioRx_2peaks_control_data 09152017.xlsx','HK2_1_062714_LBp5',2,101);
%
%   See also XLSREAD.

% Auto-generated by MATLAB on 2018/08/23 10:21:33

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 4
    endRow = 102;
end

%% Import the data
data = xlsread(workbookFile, sheetName, sprintf('A%d:J%d',startRow(1),endRow(1)));
for block=2:length(startRow)
    tmpDataBlock = xlsread(workbookFile, sheetName, sprintf('A%d:J%d',startRow(block),endRow(block)));
    data = [data;tmpDataBlock]; %#ok<AGROW>
end

%% Create table
tableout = table;

%% Allocate imported array to column variable names
tableout.No = data(:,1);
tableout.VarName2 = data(:,2);
tableout.Pyr_position = data(:,3);
tableout.Pyr_area = data(:,4);
%tableout.PyrHyd_position = data(:,5);
%tableout.PyrHyd_area = data(:,6);
%tableout.Lac_upfield = data(:,7);
%tableout.Lac_upf_area = data(:,8);
%tableout.Lac_downfield = data(:,9);
%tableout.Lac_downf_area = data(:,10);

