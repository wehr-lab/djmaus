function result=compare_versions(version1, version2)
%compare version numbers of the form 0.6.4
%returns 1 if version1 is greater, -1 if version2 is greater, 0 if equal

 % Split the version strings into parts
    parts1 = str2double(split(version1, '.'));
    parts2 = str2double(split(version2, '.'));

    % Extend shorter version with zeros if necessary
    maxLen = max(length(parts1), length(parts2));
    if length(parts1) < maxLen
        parts1(end+1:maxLen) = 0;
    elseif length(parts2) < maxLen
        parts2(end+1:maxLen) = 0;
    end

    % Compare each part from left to right
    for i = 1:maxLen
        if parts1(i) > parts2(i)
            result = 1; % version1 is greater
            return;
        elseif parts1(i) < parts2(i)
            result = -1; % version2 is greater
            return;
        end
    end

    % If all parts are equal
    result = 0; % versions are equal
end







