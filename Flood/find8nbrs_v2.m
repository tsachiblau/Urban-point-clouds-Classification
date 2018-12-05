function nbrsInds = find8nbrs_v2(idx,M,N)

[sub1,sub2] = ind2sub_v2([M N],idx);

sub1_up = sub1 - 1;
sub1_down = sub1 + 1;
sub2_left = sub2 - 1;
sub2_right = sub2 + 1;

% ul u ur
%
% l  *  r
%
% dl d dr

% u
if sub1_up>=1
    ind_up = idx - 1;
else
    ind_up = nan;
end

% d
if sub1_down<=M
    ind_down = idx + 1;
else
    ind_down = nan;
end

% l
if sub2_left>=1
    ind_left = idx - M;
else
    ind_left = nan;
end

% r
if sub2_right<=N
    ind_right = idx + M;
else
    ind_right = nan;
end

% ul
if sub1_up>=1 && sub2_left>=1
    ind_ul = idx - M - 1;
else
    ind_ul = nan;
end

% ur
if sub1_up>=1 && sub2_right<=N
    ind_ur = idx + M - 1;
else
    ind_ur = nan;
end

% dl
if sub1_down<=M && sub2_left>=1
    ind_dl = idx - M + 1;
else
    ind_dl = nan;
end

% dr
if sub1_down<=M && sub2_right<=N
    ind_dr = idx + M + 1;
else
    ind_dr = nan;
end

%
nbrsInds = [ind_up ind_down ind_left ind_right ind_ul ind_ur ind_dl ind_dr];

return