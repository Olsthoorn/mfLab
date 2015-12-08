function pauzeer(N)
if 1
   pause
else
	%Pauzeer N seconds, showing remaining time
	N=5;
	fprintf('pausing %d seconds: ',N);
	for i=N:-1:1
 		fprintf('%d',i');
 		pause(1)
	end
	fprintf('\n\n\n');
end