function [hour,minute,second] = TimeFormat(nbSeconds)
    hour = floor(nbSeconds/3600);
    nbSeconds2 = nbSeconds-hour*3600;
    minute = floor((nbSeconds-hour*3600)/60);
    second = nbSeconds2 - minute * 60;
end