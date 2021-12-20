function []=ceain(fuel,oxid,rof,pcc,acat)

delete 'cea.inp'

input = 'problem';
dlmwrite('cea.inp',input,'-append','delimiter','','newline','pc');
input = ['rocket',' fac',' ac/at=',num2str(acat),' tcest,k=',num2str(3800)];
dlmwrite('cea.inp',input,'-append','delimiter','','newline','pc');
input = ['p,bar=',num2str(pcc),','];
dlmwrite('cea.inp',input,'-append','delimiter','','newline','pc');
input = 'react';
dlmwrite('cea.inp',input,'-append','delimiter','','newline','pc');
input = ['fuel=',fuel,' wt=',num2str(1),' t,k=',num2str(50)];
dlmwrite('cea.inp',input,'-append','delimiter','','newline','pc');
input = ['oxid=',oxid,' wt=',num2str(rof)];
dlmwrite('cea.inp',input,'-append','delimiter','','newline','pc');
input = 'output';
dlmwrite('cea.inp',input,'-append','delimiter','','newline','pc');
input = ['massf',' short',' trace='];
dlmwrite('cea.inp',input,'-append','delimiter','','newline','pc');
input = ['plot',' o/f',' p',' t'];
dlmwrite('cea.inp',input,'-append','delimiter','','newline','pc');
input = ['end'];
dlmwrite('cea.inp',input,'-append','delimiter','','newline','pc');
% type cea.inp