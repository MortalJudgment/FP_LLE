function o=asknu(s,i)
% ASKNU: little function to ask user for number, with default. 
%
%  o=asks(s,i)
% s=string, i=default
as=[s,'(default=',num2str(i),'):'];
reply= str2num(input(as , 's')); 
if isempty(reply) reply=i; end 
if isnumeric(reply) o=reply; else o=i; end 
end
