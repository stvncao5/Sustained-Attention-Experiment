function list_of_valid_names_that_we_actually_want = dirP(folderName, flag_fullSearch)
    % "MATLAB can't list a directory without all the cross-OS platform madness, and it costs how much?!?"
    % This script builds on the dir() function by taking its output and "filtering" it, generating a new and improved output.
    
    % Preassigned names of folders that we want to exclude.
    excludeList = [".", "..", ".DS_Store"];
    
    % Get the full list of all the files/folders in the specified directory name (note that this is not the full path name).
    fileListing = dir(folderName);
    filenameListing = [];
    
    % Go through each item in the directory listing and append its name unless it is on the exclusion list.
    for i = 1:length(fileListing)
        
        if any(   strcmp( fileListing(i).name , excludeList )   )
            % If the filename is on the exclusion list, then do not include it (i.e. do nothing).
        else
            filenameListing = cat( 2, filenameListing, string(fileListing(i).name) );
        end
        
    end
    
    list_of_valid_names_that_we_actually_want = filenameListing;
    
end