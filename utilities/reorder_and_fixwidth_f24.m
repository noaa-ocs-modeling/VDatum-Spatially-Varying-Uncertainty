 function reorder_and_fixwidth_f24(infile, outfile, desiredOrder)
% Reorder SAL-style fort.24 blocks and rewrite node rows with fixed widths.
% Node rows format: I10, F15.6, F15.6  -> "%10d%15.6f%15.6f"
% Usage:
%   reorder_and_fixwidth_f24('fort.24','fort.24_fixed',...
%       {'K1','K2','M2','N2','O1','P1','Q1','S2'});


    if ischar(desiredOrder), desiredOrder = cellstr(desiredOrder); end
    desiredOrder = upper(desiredOrder(:)).';


    % --- read raw & keep line-ending style ---
    fid = fopen(infile,'rb'); assert(fid>=0,'Cannot open %s', infile);
    bytes = fread(fid,[1,inf],'*uint8'); fclose(fid);
    txt = char(bytes);
    eol = detect_eol(txt);


    % split into lines preserving EOL
    lines = split_keep_eol(txt);
    nL = numel(lines);


    % --- detect headers (any line containing 'SAL'), name = first token ---
    hdrIdx = []; hdrName = {};
    for i = 1:nL
        L = strip_eol(lines{i});
        if isempty(strtrim(L)), continue; end
        if contains(upper(L),'SAL')
            tok = regexp(L,'([A-Za-z0-9]+)','match','once');
            if ~isempty(tok)
                hdrIdx(end+1)   = i;          %#ok<AGROW>
                hdrName{end+1}  = upper(tok); %#ok<AGROW>
            end
        end
    end
    assert(~isempty(hdrIdx), 'No SAL headers found.');


    % --- slice blocks with SAFE indexing ---
    nB = numel(hdrIdx);
    blocks = struct('name',[],'lines',[]);
    blocks(nB).name = '';
    for k = 1:nB
        i1 = hdrIdx(k);
        if k < nB
            i2 = hdrIdx(k+1) - 1;
        else
            i2 = nL;
        end
        blocks(k).name  = hdrName{k};
        blocks(k).lines = lines(i1:i2);
    end
    currentOrder = upper({blocks.name});


    % verify all desired exist
    missing = desiredOrder(~ismember(desiredOrder, currentOrder));
    assert(isempty(missing), 'Missing block(s): %s', strjoin(missing,' '));


    % --- write in desired order; fixed-width node rows ---
    fout = fopen(outfile,'wb'); assert(fout>=0,'Cannot open %s', outfile);
    for j = 1:numel(desiredOrder)
        bi = find(strcmpi(currentOrder, desiredOrder{j}), 1);
        BL = blocks(bi).lines;


        % keep first 4 lines of the block unchanged (header + astro + 1 + name)
        headCount = min(4, numel(BL));
        for q = 1:headCount
            fwrite(fout, strip_trailing_spaces(BL{q}), 'char');
        end


        % remaining lines: node  amp  phase  -> fixed widths
        for r = (headCount+1):numel(BL)
            row = strip_eol(BL{r});
            if isempty(strtrim(row))
                fwrite(fout, BL{r}, 'char'); % blank/odd line
                continue
            end
            row2 = strrep(row,'D','E');
            tok = regexp(strtrim(row2), '\s+', 'split');
            if numel(tok) ~= 3
                fwrite(fout, BL{r}, 'char'); % non-standard, write as-is
                continue
            end
            node = str2double(tok{1});
            amp  = str2double(tok{2});
            pha  = str2double(tok{3});
            if any(isnan([node,amp,pha]))
                fwrite(fout, BL{r}, 'char');
                continue
            end
            fprintf(fout, ' %10d%15.6f%15.6f%s', round(node), amp, pha, eol);
        end
    end
    fclose(fout);
    fprintf('Reordered & fixed-width written → %s\n', outfile);
end


% ---------- helpers ----------
function e = detect_eol(t)
    if contains(t, sprintf('\r\n')), e = sprintf('\r\n'); else, e = sprintf('\n'); end
end


function L = split_keep_eol(t)
    L = {};
    i = 1; n = length(t);
    while i <= n
        j = i;
        while j <= n && t(j) ~= char(10) && t(j) ~= char(13)
            j = j+1;
        end
        if j > n
            L{end+1} = t(i:n); %#ok<AGROW>
            break
        end
        if t(j) == char(13) && j < n && t(j+1) == char(10)
            seg = t(i:j+1); j = j+2;
        else
            seg = t(i:j);   j = j+1;
        end
        L{end+1} = seg; %#ok<AGROW>
        i = j;
    end
end


function s = strip_eol(s)
    if endsWith(s, sprintf('\r\n')), s = extractBefore(s, strlength(s)-1);
    elseif endsWith(s, sprintf('\n')), s = extractBefore(s, strlength(s));
    end
end


function s = strip_trailing_spaces(s)
    if endsWith(s, sprintf('\r\n')), e = sprintf('\r\n'); core = extractBefore(s, strlength(s)-1);
    elseif endsWith(s, sprintf('\n')), e = sprintf('\n'); core = extractBefore(s, strlength(s));
    else, e = ''; core = s;
    end
    core = regexprep(core, '\s+$', '');
    s = [core e];
end


