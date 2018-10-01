%% Homework 4
%% Hakim Lachnani
%% 7 February 2017

%% Clear Workspace

clear all
close all
clc

%% Problem 5
% Halley's Comet Plot

%% Call Function

axis square
hold on
title('Orbit of Halley Comet and the Planets (in AU)');
%sun
plot(0,0,'b');
% Mercury
makeellipse(.3871,0);
% Venus
makeellipse(.7233,0);
% Earth
makeellipse(1,0);
% Mars
makeellipse(1.524,0);
% Jupiter
makeellipse(5.203,0);
% Saturn
makeellipse(9.539,0);
% Neptune
makeellipse(19.18,0);
% Uranus
makeellipse(30.07,0);
% Pluto
makeellipse(39.944,0);
% Halley
makeellipse(17.96,.967);

%% End of Homework 4 script
