function concatenate_pdfs(pdf_list, output_pdf)
    % concatenate_pdfs(pdf_list, output_pdf)
    % This function concatenates the PDFs in the input cell array `pdf_list`
    % and creates a single PDF named `output_pdf`.
    %
    % ZE note: written by ChatGPT with the following prompt:
    % "please write a matlab script that will take a cell array of pdf
    % names and write a .tex file that concatenates all the pdfs and then
    % the script should run that .tex file to create one concatenated pdf
    % with all the individual pdfs"
    
    % Check for valid input
    if isempty(pdf_list) || ~iscell(pdf_list)
        error('pdf_list must be a non-empty cell array of PDF filenames.');
    end
    
    if nargin < 2
        output_pdf = 'concatenated_output.pdf'; % Default output file name
    end
    
    % LaTeX .tex file name
    tex_filename = 'concatenate_pdfs.tex';
    
    % Open the .tex file for writing
    fid = fopen(tex_filename, 'w');
    if fid == -1
        error('Could not open .tex file for writing.');
    end
    
    % Write LaTeX preamble
    fprintf(fid, '\\documentclass{article}\n');
    fprintf(fid, '\\usepackage{pdfpages}\n');
    fprintf(fid, '\\begin{document}\n');
    
    % Add each PDF to the .tex file
    for i = 1:length(pdf_list)
        pdf_name = pdf_list{i};
        if ~isfile(pdf_name)
            fclose(fid);
            error('File "%s" does not exist.', pdf_name);
        end
        fprintf(fid, '\\includepdf[pages=-]{%s}\n', pdf_name);
    end
    
    % End the document
    fprintf(fid, '\\end{document}\n');
    
    % Close the .tex file
    fclose(fid);
    
    % Run LaTeX to generate the concatenated PDF
    cmd = sprintf('/Library/TeX/texbin/pdflatex %s', tex_filename);
    system(cmd);
    
    % Rename the output to desired name
    movefile('concatenate_pdfs.pdf', output_pdf);
    
    % Cleanup temporary files
    delete('concatenate_pdfs.aux', 'concatenate_pdfs.log', tex_filename);
    
    fprintf('Concatenated PDF created: %s\n', output_pdf);
end
