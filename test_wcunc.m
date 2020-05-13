clear;
clc;

report_filename = 'test_result.txt';
test_path = which(mfilename);
path_parts = strsplit(test_path, '/');
test_dir = strjoin(path_parts(1 : end - 1), '/');
cd(test_dir);

n_examples = size(dir('examples/'), 1) - 2;
n_tc = 2;

cnames = {'given', 'selected'};
res = cell(n_examples, 1);
warning('off', 'all');

pool = gcp;
spmd
	warning('off', 'all');
end

% print progress bar
fprintf('Tests running...\n');
fprintf('%s\n\n', repmat('.', [1, n_examples]));

% run tests
parfor kk = 1 : n_examples
	warning('off', 'all');
	ex_data = load(sprintf('./examples/example%d', kk));
	usys = ex_data.usys;
	freq = ex_data.freq;
	for tc = 1 : n_tc
		if tc == 1
			input = {usys, freq};
		elseif tc == 2
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
			er = err.message;
		end
		if pass
			obj = info.obj;
			er = [];
			stable = isempty(info.cfreq);
		else
			obj = [];
			time = [];
			stable = [];
		end
		res{kk}.(cnames{tc}).pass = pass;
		res{kk}.(cnames{tc}).obj = obj;
		res{kk}.(cnames{tc}).time = time;
		res{kk}.(cnames{tc}).er = er;
		res{kk}.(cnames{tc}).stable = stable;
	end
	fprintf('\b|\n');
end
warning on;

% write report into file
file_id = fopen(report_filename, 'w');
for kk = 1 : n_examples
	fprintf(file_id, '%d\t', kk);
	resk = res{kk};
	tcn = fieldnames(resk);
	for tc = 1 : n_tc
		resktc = resk.(tcn{tc});
		if tc > 1
			fprintf(file_id, '\t');
		end
		fprintf(file_id, '%s freqs.\t', tcn{tc});
		if resktc.pass
			fprintf(file_id, 'pass\n');
			fprintf(file_id, '\t\ttime: %d min\n\t\tobj.: %.2f%%', ceil(resktc.time / 60),...
				resktc.obj * 100);
			if ~resktc.stable
				fprintf(file_id, ' (robustly unstable)');
			end
			fprintf(file_id, '\n');			
		else
			fprintf(file_id, ['fail\t', resktc.er, '\n']);
		end
	end
end
closeflag = fclose(file_id);
if closeflag == 0
	fprintf('Report written to %s.\n', report_filename);
elseif closeflag == -1
	error('Unable to close file %s.', report_filename);
else
	error('An uncexpected even ocurred when trying to close file %s.', report_filename);
end