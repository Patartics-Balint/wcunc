clear;
clc;

test_path = which(mfilename);
path_parts = strsplit(test_path, '/');
test_dir = strjoin(path_parts(1 : end - 1), '/');
cd(test_dir);

n_examples = size(dir('examples/')) - 2;
fprintf('Tests running...\n');

delete('test_result.txt');
diary('test_result.txt');

warning off;
for kk = 1 : n_examples - 1 % the algorithm takes forever to run for example 41
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
diary('off');
warning on;