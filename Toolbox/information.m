function output = information(input)
    switch input
        case "version"
            text = "1.0 \n";
        case "license"
            fileId = fopen("license.txt");
            text = fscanf(fileId,'%c');
            fclose(fileId);
        otherwise
            error("Requested information is not available.");
    end

    if nargout > 0
        output = text;
    else
        fprintf(text);
    end
end