function chMapT = BlackrockMap()
% blackrock 1024ch H-plate omentic to Blackrock Cerebus instance channel map
% Assuming we use Cereplex-M amplifiers. Each amplifier has 128 channels
% distributed on four 32 channel banks: A(1-32)B(33-64)C(65-96)D(97-128)
% There are 32 36-pin omnetic connectors on the H-plate, Connector 1-32.
% Here is a map of the connection of the 8x Cereplex-M to 32 Omnectics:
% O01 -> 1D
% O02 -> 1A
% O03 -> 1B
% O04 -> 1C
%
% O05 -> 2D
% O06 -> 2A
% O07 -> 2B
% O08 -> 2C
%
% O09 -> 3D
% O10 -> 3A
% O11 -> 3B
% O12 -> 3C
%
% O13 -> 4D
% O14 -> 4A
% O15 -> 4B
% O16 -> 4C
%
% O17 -> 5C
% O18 -> 5B
% O19 -> 5A
% O20 -> 5D
%
% O21 -> 6C
% O22 -> 6B
% O23 -> 6A
% O24 -> 6D
%
% O25 -> 7C
% O26 -> 7B
% O27 -> 7A
% O28 -> 7D
%
% O29 -> 8C
% O30 -> 8B
% O31 -> 8A
% O32 -> 8D

% Translate this into channel map, each column is channel #, Omnetic #, channel # on Omnetic, CereplexM
% #, CereplexM bank ID, Cerebus instance #, Cerebus instance channel #

% use the most silly way to program this part for clarity of logic. The
% code does NOT need to be optimized as it runs only once.

% cereplexM 1-4
channelN = 0;
tmpBankOrder = [4,1,2,3];
for ThisCereplexM = 1:4
    for thisOmnetic = 1:4
        for thisChannel = 1:32
            channelN = channelN + 1;
            tmpOmnetic = (ThisCereplexM-1) * 4 + thisOmnetic;
            chMap(channelN,:) = [channelN,tmpOmnetic,thisChannel,ThisCereplexM,tmpBankOrder(thisOmnetic),ThisCereplexM,(tmpBankOrder(thisOmnetic)-1)*32+thisChannel];
        end
    end
end

% cereplexM 5-8
tmpBankOrder = [3,2,1,4];
for ThisCereplexM = 5:8
    for thisOmnetic = 1:4
        for thisChannel = 1:32
            channelN = channelN + 1;
            tmpOmnetic = (ThisCereplexM-1) * 4 + thisOmnetic;
            chMap(channelN,:) = [channelN,tmpOmnetic,thisChannel,ThisCereplexM,tmpBankOrder(thisOmnetic),ThisCereplexM,(tmpBankOrder(thisOmnetic)-1)*32+thisChannel];
        end
    end
end

% make it into a table
chMapT.channelNumber = chMap(:,1);
chMapT.OmneticConnectorNumber = chMap(:,2);
chMapT.channelOnOmnetic = chMap(:,3);
chMapT.CereplexMnumber = chMap(:,4);
chMapT.CereplexMbankID = chMap(:,5);
chMapT.CerebusInstanceNumber = chMap(:,6);
chMapT.CerebusInstanceChannelNumber = chMap(:,7);