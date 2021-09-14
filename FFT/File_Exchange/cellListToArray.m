function [array] = cellListToArray(cellList,columnNumber)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%arraySize = [max(cellfun('length',cellList(:,columnNumber))),size(cellList,1)];
array = cell2mat(cellList(1,columnNumber));
for i =2:size(cellList,1)
    array = padconcatenation(array,cell2mat(cellList(i,columnNumber)),2);
end

