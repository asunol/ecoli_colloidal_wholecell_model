%generates random numbers in interval using uniform distribution
%in: 
%   [a,b] limits
%   M:  cols
%   N:  rows
%out: MxN array 
    
function r = randInterval (a,b,M,N)
 r = a + (b-a).*rand(M,N);
end





