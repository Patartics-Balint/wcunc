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

for kk = 1 : n_examples
	warning('off', 'all');
	fprintf('%d\t', kk);
	load(sprintf('./examples/example%d', kk));
	for tc = 1 : 2
		if tc == 1
			fprintf('given freqs.\t');
			input = {usys, freq};
		elseif tc == 2
			fprintf('\tselected freqs.\t');
			input = {usys};
		else
			error('Test case not recognised.');
		end
		try
			tic;
			[wcu, wcsys, info] = wcunc(input{:});
			time = toc;
			pass = true;
		catch err
			pass = false;
		end
		if pass
			fprintf('pass\n');
			obj_th = wcgainlbgrid(usys, info.freq);
			obj_th = sum(obj_th);
			obj = sigma(wcsys, info.freq);
			obj = sum(obj(1, :));
			fprintf('\t\ttime: %d sec\n\t\tobj.: %.2f%%\n', ceil(time), obj / obj_th * 100);
		else
			fprintf(['fail\t', err.message, '\n']);
		end
	end
end
diary('off');
warning on;

% remove warnings form report
report = fileread('test_result.txt');
no_backspaces = regexp(report, '\b', 'split');
report = strjoin(no_backspaces, '');
no_warnings = regexp(report, '\[[^\]]+\]', 'split');
report = strjoin(no_warnings, '');
no_empty_lines = regexp(report, ' \n', 'split');
report = strjoin(no_empty_lines, '');

% write report into file
file_id = fopen('test_result.txt', 'w');
fprintf(file_id, report);