function data = read_edi_csv(filename)
%READ_EDI_CSV Read an EDI answers-by-gender CSV file.
%
%   DATA = READ_EDI_CSV(FILENAME)
%
%   Reads a CSV file with the following structure:
%
%       # Section
%       # Question
%       option,female,male
%       Never,21,45
%       Rarely/single instance,3,2
%       ...
%
%   The returned structure contains:
%       data.filename   - input filename
%       data.section    - section/title from the first comment line
%       data.question   - question from the second comment line
%       data.answers    - table containing the actual data
%
%   Example:
%
%       data = read_edi_csv('edi-answers-by-gender-3-03.csv');
%
%       data.answers
%       data.question

    % Open file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    cleanup = onCleanup(@() fclose(fid));

    % Read first two comment lines
    line1 = fgetl(fid);
    line2 = fgetl(fid);

    if ~ischar(line1) || ~ischar(line2)
        error('File does not contain the expected header information.');
    end

    % Remove "# " from the beginning
    section = strtrim(erase(line1, '#'));
    question = strtrim(erase(line2, '#'));

    % Read the remaining CSV data
    % The first two lines have already been consumed.
    opts = detectImportOptions(filename, ...
        'Delimiter', ',', ...
        'NumHeaderLines', 2);

    answers = readtable(filename, opts);

    % Store everything in a structure
    data = struct();
    data.filename = filename;
    data.section  = section;
    data.question = question;
    data.answers  = answers;

end