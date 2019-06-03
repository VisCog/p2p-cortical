% Winawer_datasummary

flist = [];
durlist = [];
pwlist = [];
for i=1:5
    load(['trial_data_site', num2str(i), '.mat'])
    flist = cat(2, flist, unique([frequency.val])); flist = unique(flist);
    durlist = cat(2, durlist, unique([duration.val])); durlist = unique(durlist);
    pwlist = cat(2, pwlist, unique([pulsewidth.val])); pwlist = unique(pwlist);
end
mat = zeros(length(flist), length(durlist), length(pwlist));
for i=1:5
for ff = 1:length(flist)
        for dd = 1:length(durlist)
            for pp = 1:length(pwlist)
                ind = find([[frequency.val] == flist(ff)] & [[duration.val] == durlist(dd)] & [[pulsewidth.val] ==pwlist(pp)])
                mat(ff, dd, pp) = mat(ff, dd, pp) + length(ind);
            end
        end
end
end


       

    