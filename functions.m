%function [outputArg1,outputArg2] = untitled(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%outputArg1 = inputArg1;
%outputArg2 = inputArg2;
%end

function pounds = KG2P(kg)
    pounds = kg / 0.453592;
end
function kg = P2KG(pounds)
    kg = pounds * 0.453592;
end
function inch = CM2INCH(cm)
    inch = cm/2.54;
end
function cm = INCH2CM(inch)
    cm = inch*2.54;
end

