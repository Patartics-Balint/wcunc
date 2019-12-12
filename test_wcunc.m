clear;
clc;
n_examples = size(dir('examples/')) - 2;
fprintf('Tests running...\n');
warning off;
for kk = 1 : n_examples
	fprintf('%d\t', kk);
	load(sprintf('./examples/example%d', kk));
	fprintf('given freqs.\t');
	try 
		[wcu, wcg, info] = wcunc(usys, freq);
		fprintf('pass\n');
	catch err
		fprintf(['fail\t', err.message, '\n']);
	end
	fprintf('\tselected freqs.\t');
	try 
		[wcu, wcg, info] = wcunc(usys);
		fprintf('pass\n');
	catch err
		fprintf(['fail\t', err.message, '\n']);
	end	
end
warning on;