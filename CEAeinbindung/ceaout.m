function [C]=ceaout

fid = fopen('cea.out');

C = textscan(fid,'%s');