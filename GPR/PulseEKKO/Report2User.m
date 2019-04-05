function Report2User(reportMessage,filename,workingDirectory)
	
	% Open the Report Log file for appending.
	fullFileName = fullfile(workingDirectory,filename);
	fid = fopen(fullFileName, 'at');
	fprintf(fid, '%s\n', reportMessage); % To file
	fclose(fid);
	return; % from Report2User()
end