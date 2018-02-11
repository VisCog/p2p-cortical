function ebs_downloadData()
% Download and unzip the data from the Open Science Framework project page
% associated with this paper:
%
%   Winawer and Parvizi (2016). Linking Electrical Stimulation of Human
%   Primary Visual Cortex, Size of Affected Cortical Area, Neuronal
%   Responses, and Subjective Experience Neuron. 92(6):1213?1219
%   http://dx.doi.org/10.1016/j.neuron.2016.11.008
%
% Alternatively, the data can be downloaded manually from this site:
% https://osf.io/zgvvm/
% Or from this DOI: 10.17605/OSF.IO/PZ42U
%
% The code downloads a single zip file (41.3 MB), places it in the root
% directory of the project, and unzips it into the folder named 'data'

url = 'https://osf.io/zgvvm/?action=download';
pth = fullfile(ebsRootPath, 'ebs_data.zip');
fname = websave(pth, url);
unzip(fname);

end