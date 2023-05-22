function [id] = HGSid(species)
% HGS - Thermochemistry of gas mixture
%
% Property of THRUST, unauthorized distribution is not allowed
% version: HGSid 2.2
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % Branched and modified from the original ManelSoria/HGS repository %
% % Original head:                                                    %
% % *HGS 2.1                                                          % 
% % *By Caleb Fuster, Manel Soria and Arnau Miró                      %
% % *ESEIAAT UPC                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%**************************************************************************
%
% [id] = HGSid(species)
%
%**************************************************************************
% 
% HGSid finds the id of the species to improve the speed of the code
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or numbers of species
%
% Outputs:
%--------------------------------------------------------------------------
% id --> Species id
%
%**************************************************************************

%{
Changelog:
  > version: 2.2 - 14/12/2022 - Alessandro Rampazzo
    - moved the "error('uhh ? HGSid wrong data type')" up the nested
    function find1
    - now it accepts also strings and string arrays as species

  > version: 2.1 - 14/12/2022 - Alessandro Rampazzo
     - branched from the original ManelSoria/HGS repository
%}


if isnumeric(species)
    id=species;
    return
end

global HGSdata;HGSload;


ns=size(HGSdata.Mm,1); % Number of species
if isfield(HGSdata,'comb') && ~isempty(HGSdata.comb)
    combination = 1;
else
    combination = 0;
end

if iscell(species)    
    id = zeros(1,length(species));
    for ii = 1:length(species)
        id(ii) = find1(species{ii});
    end
    return

elseif ischar(species)
    id = find1(species);
    return

elseif isstring(species)
    if length(species) == 1
        id = find1(char(species));
        return
    else
        id = zeros(1,length(species));
        for ii = 1:length(species)
            id(ii) = find1(species(ii));
        end
    end
    return
end

error('uhh ? HGSid wrong data type')

function id = find1(name)
    search = ismember(HGSdata.name(:),name);
    if combination
        search1 = ismember(HGSdata.comb(:),name);
    else
        search1 = 0;
    end

    if all(search == 0) && all(search1 == 0)
        error('HGSid:  %s not found in the data base\n',name);
    elseif ~all(search == 0)
        id = find(search);
    else
        id = ns + find(search1);
    end
end

end