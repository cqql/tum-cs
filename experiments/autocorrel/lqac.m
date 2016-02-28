function [p s r b] = main(f)
	t1 = cputime;
	p = loadlq(f);
	t2 = cputime;
	disp(['loaded in ' num2str(t2-t1)]);
	[dy dx] = size(p);
	b = blowup(p(:,3));
	s = tdiffs(p, 1, 3);
	r = tdiffs(p, 2, 3);
	t3 = cputime;
	disp(['processed in ' num2str(t3-t2)]);
	cdfplot(s)
	print([f 'ps'],'-dpng')
	cdfplot(r)
	print([f 'pr'],'-dpng')
	autocorr(b,round(dy/10));
	print([f 'pc'],'-dpng')
	disp(['plotted in ' num2str(cputime-t3)]);
end

function d = tdiffs(p,sc,dc)
	% I bet matlab has a builtin for this…
	[r c] = size(p);
	ret = [];
	for i = 2:r
		td = double(p(i,sc) - p(i-1,sc)) / 1e9;
		sd = double(p(i,dc) - p(i-1,dc));
		ret = [ret; td / sd];
	end
	d = ret;
end

function r = loadlq(path)
	t1 = cputime;
	read = cell(0,1);
	fid = fopen(path);
	tline = fgets(fid);
	while ischar(tline)
		read{end+1,1} = tline;
		tline = fgets(fid);
	end
	fclose(fid);
	t2 = cputime;
	disp(['read in ' num2str(t2-t1)]);
	l = size(read,1);
	res1 = cell(l,1);
	res2 = cell(l,1);
	res3 = cell(l,1);
	t3 = cputime;
	disp(['allocated in ' num2str(t3-t2)]);
	parfor i = 1:l
		s = sscanf(read{i}, '%u-%u %u-%u %s %u');
		if size(s,1) == 6
			res1{i,1} = totime(s([1 2]));
			res2{i,1} = totime(s([3 4]));
			res3{i,1} = int64(s(6));
		else
			disp(['Warning, ignoring line ' num2str(i) ': "' strtrim(tline) '"']);
		end
	end
	t4 = cputime;
	disp(['parsed in ' num2str(t4-t3)]);
	r = [cell2mat(res3) cell2mat(res2) cell2mat(res3)];
	disp(['reallocated in ' num2str(cputime-t4)]);
end

function t = totime(s)
	fs = s(1);
	ns = s(2);
	t = int64(fs) * int64(1e9) + int64(ns);
end

function b = blowup(p)
	p = [0; p];
	[r c] = size(p);
	ret = [];
	for i = 2:r
		fl = ones(p(i) - p(i-1) - 1,1) * -1;
		ret = [ret; fl; 1];
	end
	b = ret;
end
	
