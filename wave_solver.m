%This matlab code can be used as verification for the C MPI program
%the code is taken from http://cs.slu.edu/~chambers/fall09/cs145/wave-equation.pdf

format long
dx = .00100;
dt =.0001;
c = 4;
L = 6;
stopTime = 5;
r = c*dt/dx;
n = L/dx
current = .5-.5*cos(2*pi/L*[0:dx:L]);
past = current;

for t=0:dt:stopTime
    future(1)=0;
    future(2:n-1)=r^2*(current(1:n-2)+current(3:n))+2*(1-r^2)*current(2:n-1)-past(2:n-1);
    future(n)=0;
    past = current;
    current = future;
%         
%     if mod(t/dt,100)==0
%         plot(current)
%         axis([0 6001 -1 1])
%         future(n/2)
%         pause(.001)
%     end
end

plot(current)
axis([0 6001 -1 1])

    
