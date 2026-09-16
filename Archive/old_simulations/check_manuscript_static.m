function nfail = check_manuscript_static(fname)
% nfail = check_manuscript_static(['notebooks/smi_manuscript_60deg.m'])
%
% Static checks on the manuscript simulation, so a structural mistake is caught
% in seconds rather than after a multi-hour run.
%
%   octave-cli --no-gui -q check_manuscript_static.m
%   >> check_manuscript_static
%
% It does NOT run the simulation. It answers the narrower question "will this
% file get to the end without throwing", which is the question worth answering
% before handing a run to somebody's overnight machine.
%
% WHAT IT CHECKS, AND WHY EACH ONE EXISTS
%
% 1. *It parses.* The whole file is wrapped in |if false ... end| and eval'd, so
%    every line is parsed and none is executed. Catches syntax errors anywhere,
%    including inside branches a short run would never enter.
%
% 2. *Subscript arity on the scoring arrays.* This is the one that earned its
%    place. When the CSD arms were wired in, |res_all| and friends gained a
%    leading ARM index and became 4-D. Every read had to gain a subscript, and
%    one did not -- the final summary table, 500 lines below the change, which
%    only executes after everything else has already printed. A 12-minute run
%    ended in `out of bound` there. This check reads every subscripted use of
%    those arrays and reports any with the wrong count.
%
%    The comma counting is BRACKET-DEPTH AWARE, which is not optional: a naive
%    version counts the comma inside |max(a, [], 2)| or |sums{im}(ia, iL, ...)|
%    and reports nonsense. "README for Claude" section 4 records that a naive
%    keyword balance check is useless for exactly this reason.
%
% 3. *Every RUN{} field that is read is also written.* The per-kernel results
%    are stashed in a struct built in one place and read by the figures. A typo
%    or a renamed field there fails only when the figures run, at the very end.
%
% Exit status is the number of failures, so it can gate a commit.

if nargin < 1 || isempty(fname)
    here = fileparts(mfilename('fullpath'));
    if isempty(here), here = pwd; end
    fname = fullfile(here, 'notebooks', 'smi_manuscript_60deg.m');
end
if ~exist(fname, 'file'), error('no such file: %s', fname); end

VERDICT = {'** FAILED **', 'ok'};
nfail = 0;
src   = fileread(fname);
lines = strsplit(src, sprintf('\n'));

fprintf('\n=== static checks on %s ===\n', fname);

% ---------------------------------------------------------------- 1. parse
ok_parse = true; msg = '';
try
    eval(['if false,' sprintf('\n') src sprintf('\n') 'end']);
catch err
    ok_parse = false; msg = err.message;
end
fprintf('   CHECK the file parses                                  %s\n', ...
        VERDICT{1+ok_parse});
if ~ok_parse
    fprintf(2, '        %s\n', msg);
    nfail = nfail + 1;
    fprintf('\n   Stopping: nothing below is meaningful if it does not parse.\n');
    return
end

% ------------------------------------------------- 2. scoring array arity
% name -> how many subscripts a read of it must have.
arrays = { 'res_all',  4; 'bias_all', 4; 'sd_all',   4; ...
           'med_all',  4; 'spur_all', 4; 'amp_all',  4; 'iso_all',  4; ...
           'ceil_n',   3; 'ceil_err', 3 };

% Aliases matter, and missing them is how this check first failed its own test.
% The bug that motivated it reads the arrays through a CELL built from them --
%
%     sums = {res_all, bias_all, sd_all, spur_all};
%     ... sums{im}(iL, SNR_ORD(k), ic)          % three subscripts, needs four
%
% -- so a checker that only looks for `res_all(` sails straight past it. Any
% cell literal whose members are all tracked arrays of the same arity is
% registered here as an alias with that arity, and `alias{..}(..)` is then
% checked like a direct read.
for i = 1:numel(lines)
    ln  = strip_comment(lines{i});
    tok = regexp(ln, '^\s*([A-Za-z_]\w*)\s*=\s*\{([^{}]*)\}\s*;', 'tokens', 'once');
    if isempty(tok), continue, end
    members = strtrim(strsplit(tok{2}, ','));
    ar = [];
    for m = 1:numel(members)
        j = find(strcmp(arrays(:,1), members{m}), 1);
        if isempty(j), ar = []; break, end
        ar(end+1) = arrays{j,2}; %#ok<AGROW>
    end
    if ~isempty(ar) && all(ar == ar(1))
        arrays(end+1,:) = {tok{1}, ar(1)}; %#ok<AGROW>
        fprintf('   alias: %s = {%s} tracked at %d subscripts\n', ...
                tok{1}, tok{2}, ar(1));
    end
end

fprintf('   subscript arity on the arm-indexed scoring arrays:\n');
for a = 1:size(arrays,1)
    nm   = arrays{a,1};
    want = arrays{a,2};
    bad  = {};
    for i = 1:numel(lines)
        ln = strip_comment(lines{i});
        if isempty(strfind(ln, nm)), continue, end
        uses = subscript_uses(ln, nm);
        for u = 1:numel(uses)
            got = uses{u}.n;
            if got ~= want
                bad{end+1} = sprintf('line %d: %s%s(%s)  -- %d subscripts, expected %d', ...
                                     i, nm, uses{u}.cell, uses{u}.txt, got, want); %#ok<AGROW>
            end
        end
    end
    fprintf('     %-10s expects %d   %s\n', nm, want, VERDICT{1+isempty(bad)});
    for k = 1:numel(bad)
        fprintf(2, '        %s\n', bad{k});
    end
    if ~isempty(bad), nfail = nfail + 1; end
end

% ------------------------------------------- 3. RUN{} fields read vs written
% Written as struct('name', ..., 'K', ..., ...) in one place; read as
% RUN{ikk}.field all over the figures.
written = {};
tok = regexp(src, 'RUN\{ik\}\s*=\s*struct\((.*?)\);', 'tokens', 'once');
if ~isempty(tok)
    fields = regexp(tok{1}, '''([A-Za-z_]\w*)''\s*,', 'tokens');
    for k = 1:numel(fields), written{end+1} = fields{k}{1}; end %#ok<AGROW>
end
read = {};
r = regexp(src, 'RUN\{[^}]*\}\.([A-Za-z_]\w*)', 'tokens');
for k = 1:numel(r), read{end+1} = r{k}{1}; end %#ok<AGROW>
read = unique(read);
missing = setdiff(read, written);
fprintf('   CHECK every RUN{} field read is also written           %s\n', ...
        VERDICT{1+isempty(missing)});
if ~isempty(missing)
    for k = 1:numel(missing)
        fprintf(2, '        read but never written: RUN{}.%s\n', missing{k});
    end
    nfail = nfail + 1;
else
    fprintf('        %d fields written, %d read\n', numel(written), numel(read));
end

fprintf('\n');
if nfail == 0
    fprintf('All static checks passed. This does NOT mean the numbers are right --\n');
    fprintf('it means the file should reach the end. Run it to get results.\n\n');
else
    fprintf('%d static check(s) FAILED.\n\n', nfail);
end
end

% =====================================================================
function s = strip_comment(s)
% Drop a trailing % comment, ignoring % inside single-quoted strings. Crude but
% adequate: this file has no transpose-vs-quote ambiguity on comment lines.
inq = false;
for i = 1:numel(s)
    if s(i) == ''''
        inq = ~inq;
    elseif s(i) == '%' && ~inq
        s = s(1:i-1); return
    end
end
end

% =====================================================================
function uses = subscript_uses(ln, nm)
% Every `nm(...)` and `nm{...}(...)` in ln, with the top-level comma count of
% the PAREN subscript list.
%
% Depth aware: commas inside nested (), [] or {} do not count, so
% `sums{im}(ia, iL, SNR_ORD(k), ic)` reads as 4 and not 5. The `{...}` form is
% handled because the arrays are also read through a cell alias, which is
% precisely where the bug this file exists to catch was hiding.
uses = {};
idx  = strfind(ln, nm);
for t = 1:numel(idx)
    p = idx(t);
    % must be a whole identifier, not the tail of a longer one
    if p > 1 && (isletter(ln(p-1)) || ln(p-1) == '_' || any(ln(p-1) == '0123456789'))
        continue
    end
    q = p + numel(nm);
    if q > numel(ln), continue, end
    % an identifier character straight after means this is a longer name
    if isletter(ln(q)) || ln(q) == '_' || any(ln(q) == '0123456789'), continue, end
    cellpart = '';
    if ln(q) == '{'                    % cell index first: nm{..}(..)
        d = 0; r = q;
        while r <= numel(ln)
            if any(ln(r) == '([{'), d = d + 1;
            elseif any(ln(r) == ')]}'), d = d - 1; if d == 0, break, end
            end
            r = r + 1;
        end
        if r > numel(ln), continue, end
        cellpart = ln(q:r);
        q = r + 1;
    end
    if q > numel(ln) || ln(q) ~= '(', continue, end
    open = q;                          % position of '('
    depth = 0; q = open; close = -1;
    while q <= numel(ln)
        c = ln(q);
        if any(c == '([{')
            depth = depth + 1;
        elseif any(c == ')]}')
            depth = depth - 1;
            if depth == 0, close = q; break, end
        end
        q = q + 1;
    end
    if close < 0, continue, end        % unterminated on this line; skip
    inner = ln(open+1:close-1);
    if isempty(strtrim(inner)), continue, end
    depth = 0; ncomma = 0;
    for q = 1:numel(inner)
        c = inner(q);
        if any(c == '([{'), depth = depth + 1;
        elseif any(c == ')]}'), depth = depth - 1;
        elseif c == ',' && depth == 0, ncomma = ncomma + 1;
        end
    end
    uses{end+1} = struct('n', ncomma + 1, 'txt', strtrim(inner), ...
                         'cell', cellpart); %#ok<AGROW>
end
end
