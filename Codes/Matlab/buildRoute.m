function route = buildRoute2(routeLinkIds, networkId, startTime, ...
                             endTime, direction, rho_max, startHour, ...
                             endHour)
    
    
    javaaddpath(fullfile('lib', 'CORE-release.jar'));
    javaaddpath(fullfile('lib', 'NETCONFIG.jar'));
    addpath('data_functions');

    import netconfig.*;
    import core.*;

    % Errors:
    % 1. sensor offsets (this might result in sensors being placed
    % in the wrong cells - maybe one or two more than the expected
    % cell - not a big deal)
    % 2. "hours, minutes, seconds" when constructing the
    % lineToSelect list. The hours considered are the first 14
    % hours of the day, not from 7AM to 8PM as expected (not a big
    % deal either)
    
    if nargin == 0    
        %I880 northbound
        routeLinkIds  = [186836 315387 1238 150305 233028 97544 47187 129727 191223 ...
                    191222 47190 259805 97540 190831 190830 41558 320521 320522 320523 ...
                    150310 150311 129730 83134 151055 101004 278401 150937 150944 155637 ...
                    259138];%ids of model graph links
        networkId = 132;
        
        startTime = '2012-03-05 00:00:00.000';
        endTime = '2012-03-06 00:00:00.000';
        startHour = 7;
        endHour = 20;
        direction = 'N';
        rho_max = 1/7;
    end
    
    core.Monitor.set_nid(networkId);
    net = netconfig.Network();

    % Retrieves the link info (in no particular order)
    linkInfo = retrieveLinkInfo(routeLinkIds, net);
    
    % Retrieves the sensor info (in no particular order)
    sensorInfo = retrieveSensorInfo(routeLinkIds, net, direction);
    
    % Retrieves the filtered data (ordered by date and pems id)
    filteredData = retrieveFilteredData(...
        sensorInfo(1).pemsIds, ...
        startTime, endTime, startHour, endHour);
    
    % orders the link ids according to routeLinkIds, reduces the
    % sensor info struct to only contain live sensors (i.e. only sensors
    % that show activity in the date range considered) , and orders
    % the filtered data by date and according to the order of
    % routeLinkIds
    
    linkInfo = orderStruct(linkInfo, routeLinkIds);
    
    sensorInfo = ...
        orderStruct(...
            shrinkStruct(sensorInfo(1), ...
                         unique(filteredData(1).pemsIds)), ...
            routeLinkIds);
    
    [cellLength, cellPosition, observationMatrix, sensorCellMap] = ...
        constructRouteFromLinks(linkInfo, sensorInfo);
    

    dx_min = min(cellLength);
    
    %careful, it is the assimilation time step, i.e the time 
    %between measurements. We are using filtered PEMS data
    % so dtAssimilation is 30 seconds
    dtAssimilation = 30;
    
    nbCells = sum(linkInfo.numCells);
    nbSensors = length(sensorInfo(1).pemsIds);
    
    nbHours = endHour - startHour + 1;

    totalSec = 3600*nbHours-30;%one day
    nbEnkfSteps = totalSec/30 + 1;
    sec = (startHour*3600) : 30 : (startHour*3600 + totalSec);

    %H is not the same at each time step
    activeSensors = cell(nbEnkfSteps, 1);
    densityMeasured = zeros(nbCells,nbEnkfSteps);

    for iTime = 1:nbEnkfSteps
        [hour,minute,second] = secondsToHoursMinutesSeconds(sec(iTime));
        indexes1 = find(filteredData(1).hour == hour);
        indexes2 = find(filteredData.minute(indexes1) == minute);
        indexes3 = find(filteredData.second(indexes1(indexes2)) == second);
        indexes = indexes1(indexes2(indexes3));
        activePemsIds = filteredData.pemsIds(indexes);
        numActive = size(activePemsIds,1);
        selectSensors = false(1,nbSensors);
        
        for iActive = 1:numActive
            indexSensor = find(sensorInfo(1).pemsIds == ...
                               activePemsIds(iActive));
            selectSensors(indexSensor) = true;
            densityMeasured(sensorCellMap(indexSensor), iTime) = ...
                min(filteredData.density(indexes(iActive)), ...
                    rho_max);
        end
        activeSensors{iTime} = find(selectSensors == true)';
    end

    closeNetwork(net);
    
    % OUTPUT:
    route = struct();
    
    route.filteredData = filteredData;
    route.sensorInfo = sensorInfo;
    route.linkInfo = linkInfo;    
    route.startHour = startHour;
    route.endHour = endHour;
    route.cellLength = cellLength;
    route.cellPosition = cellPosition;
    route.totalSec = totalSec;
    route.nbCells = nbCells;
    route.activeSensors = activeSensors;
    route.observationMatrix = observationMatrix;
    route.densityMeasured = densityMeasured;
    route.sensorCellMap = sensorCellMap;
    route.dx_min = dx_min;
 end


%% constructs cells from model graph links and maps sensor
%% locations to cells
function [cellLength, cellPosition, observationMatrix, sensorCellMap] = ...
        constructRouteFromLinks(linkInfo, sensorInfo)

    nbSensors = length(sensorInfo.pemsIds);
    nbLinks = length(linkInfo.linkIds);
    nbCells = sum(linkInfo.numCells);
    linkStart = zeros(1,nbLinks);
    cellLength=zeros(1,nbCells);
    cellPosition=zeros(1,nbCells);
    
    counter = 1;
    for i = 1:nbLinks
        linkStart(i)= counter;
        cellLength(1,counter:counter + linkInfo.numCells(i) - 1)= ...
            linkInfo.linkLengths(i) / linkInfo.numCells(i);
        for j = counter : (counter + linkInfo.numCells(i) - 1)
            %position is at the beginning of the cell
            %maybe it is better to put at the center of the cell?
            if j > 1
                cellPosition(j) = cellPosition(j-1) + cellLength(j-1);
            end
        end
        counter = counter + linkInfo.numCells(i);
    end
    
    sensorCellMap = zeros(1,nbSensors);
    linkSensorMap = mapMGLinksToSensors(linkInfo, sensorInfo);
    
    for i=1:nbSensors
        sensorCellMap(i) = linkStart(linkSensorMap(i))+ ...
            floor(sensorInfo.offset(i)/ ...
                  linkInfo.linkLengths(linkSensorMap(i)) * ...
                  linkInfo.numCells(linkSensorMap(i)));
    end
    
    observationMatrix = zeros(nbSensors,nbCells);
    
    for s=1:nbSensors
        observationMatrix(s,sensorCellMap(s)) = 1;
    end 
    
    %constant observation matrix (does not depend on time step n)
    %number of rows=nbSensors, number of columns=nbCells 
    %encodes the location of the sensors : H(i,j) = 1 if sensor i is in cell j,
    %0 otherwise
    %TODO:look if H is sparse: in this case should we use sparse Matlab tools
    %(look-ahead for when we want to implement it on networks)
end

function closeNetwork(net)
    net.close();
    clear net;
end

%% does as advertised :P
function linkSensorMap = mapMGLinksToSensors(linkInfo, sensorInfo)
    nbSensors = length(sensorInfo.pemsIds);
    linkSensorMap = zeros(1, nbSensors);
    
    for i=1:nbSensors
        linkSensorMap(i) = find(linkInfo.linkIds==sensorInfo.linkIds(i));
    end
end

%% retrieves the link info
function linkData = retrieveLinkInfo(routeLinkIds, net)
    
    linkPSName = 'get link info using link id list';
    linkPS = ...
        ['SELECT DISTINCT ON (nl.fk_lid) ' ...
         , ' nl.fk_lid, nl.nb_cells ' ...
         , ' FROM model_graph.network_links nl ' ...
         , ' JOIN model_graph.link_navteq_link lnl ' ...
         , ' ON (nl.fk_lid = lnl.fk_lid) ' ...
         , ' WHERE nl.fk_lid = ANY(?) AND nl.fk_nid = ?;'];
    
    linkColumnTypes = {'int', 'short'};
    jLinkIds = javaArray('java.lang.Integer', length(routeLinkIds));
        
    for iLink = 1:length(routeLinkIds)
        jLinkIds(iLink) = java.lang.Integer(routeLinkIds(iLink));
    end
        
    linkSetData = {jLinkIds, java.lang.Integer(net.nid)};
    linkSetTypes = {'aint', 'int'};

    
    linkFields = {'linkIds', 'numCells'};
    linkData = ...
        convertToMatlabStruct(...
            collectDataFromDB(...
                linkPSName, linkPS, linkColumnTypes, linkSetData, linkSetTypes), ...
            linkFields);
    
    
    linkLengths = calculateMGLinkLength(linkData.linkIds, net);
    linkData = setfield(linkData, 'linkLengths', linkLengths);
end


%% retrieves the sensor information
function sensorInfo = retrieveSensorInfo(linkIds, net, direction)
    
    sensorPSName = 'get sensor data using link id list';
    sensorPS = ...
        ['SELECT DISTINCT ON (p.id) ' ...
         , ' p.id, p.off_set, p.fk_link_id, lnl.fk_lid, lnl.ndx ' ... 
         , ' FROM pems.prop p ' ...
         , ' JOIN model_graph.link_navteq_link lnl ' ...
         , ' ON (p.fk_link_id = lnl.fk_link_id) ' ...
         , ' JOIN pems.vds_config v ' ...
         , ' ON (v.id = p.fk_id) ' ...
         , ' WHERE lnl.fk_lid = ANY(?) AND v.freeway_dir = ?; '];
    
    sensorColumnTypes = {'int', 'real', 'long', 'int', 'int'};
    jLinkIds = javaArray('java.lang.Integer', length(linkIds));
        
    for iLink = 1:length(linkIds)
        jLinkIds(iLink) = java.lang.Integer(linkIds(iLink));
    end
    
    sensorSetData = {jLinkIds, direction};
    sensorSetTypes = {'aint', 'string'};
    
    sensorFields = {'pemsIds', 'offset', 'navteqLinkIds', ...
                    'MGLinkIds', 'index'};

    sensorData = ...
        collectDataFromDB(...
            sensorPSName, sensorPS, sensorColumnTypes, ...
            sensorSetData, sensorSetTypes);
    
    sensorInfo = calculateMGLinkOffset(sensorData, net);
end

%% converts a matlab cell to a matlab structure
function mStruct = convertToMatlabStruct(matlabCell, structFields)
    [numAttributes, numFields] = size(matlabCell);
    
    mStruct = struct();
    
    for iField = 1:numFields
        if strcmp(structFields{iField}, 'date')
            attributeVector = zeros(numAttributes, 6);
            
            for iAttribute = 1:numAttributes
                [attributeVector(iAttribute, 1), ...
                 attributeVector(iAttribute, 2), ...
                 attributeVector(iAttribute, 3), ...
                 attributeVector(iAttribute, 4), ...
                 attributeVector(iAttribute, 5), ...
                 attributeVector(iAttribute, 6)] = ...
                    javaTimeToHoursMinutesSeconds(...
                        matlabCell{iAttribute, iField});
            end
            
            timeFields = {'year', 'month', 'day', ...
                          'hour', 'minute', 'second'};
            
            for iTimeField = 1:length(timeFields)
                mStruct = ...
                    setfield(...
                        mStruct, timeFields{iTimeField}, ...
                        attributeVector(:, iTimeField));
            end
            
        else
            attributeVector = zeros(numAttributes, 1);
            for iAttribute = 1:numAttributes
                attributeVector(iAttribute) = ...
                    matlabCell{iAttribute, iField};
            end
            
            mStruct = setfield(mStruct, structFields{iField}, ...
                                        attributeVector);
        end
        
    end
    
end

%% calculates the length of model graph links in the network with
%% ID nid
function linkLengths = calculateMGLinkLength(linkIds, net)
    numLinks = length(linkIds);
    linkLengths = zeros(numLinks, 1);

    for iLink = 1:numLinks
        linkLengths(iLink) = ...
            net.getLinkWithID(linkIds(iLink)).getLength();
    end
end

function sensorInfo = calculateMGLinkOffset(sensorData, net)
    numLinks = length(sensorData);
    sensorInfo = struct();
    
    for iLink = 1:numLinks
        navteqLinks = ...
            net.getLinkWithID(sensorData{iLink, 4}.intValue()) ...
            .getNavteqLinks();
        
        navLink = navteqLinks(sensorData{iLink, 5}.intValue() + 1);
        
        navSpot = ...
            netconfig.Spot(...
                navLink, sensorData{iLink, 2}.floatValue(), 0);
        
        mgSpot = navSpot.toModelGraphSpot(net);
        sensorInfo.pemsIds(iLink) = sensorData{iLink, 1};
        sensorInfo.linkIds(iLink) = mgSpot.link.id;
        sensorInfo.offset(iLink) = mgSpot.offset;
    end
    
end

%% get the filtered data for the sensor IDs passed between startTime and endTime
function filteredData = retrieveFilteredData(sensorIdList, startTime, ...
                                             endTime, startHour, endHour)
    
    
    startHour = java.lang.Integer(startHour);
    endHour = java.lang.Integer(endHour);
    
    startTime = core.Time.newTimeFromStringLikeBerkeleyDB(startTime);
    endTime = core.Time.newTimeFromStringLikeBerkeleyDB(endTime);
        
    % Sensor data over time
    filteredPSName = 'retrieve PeMS filtered data';
    filteredPS = ...
        ['SELECT f.fk_id, f.density, f.fraction_working, f.date' ...
         ' FROM pems.filtered_30sec f ' ...
         ' WHERE f.fk_id = ANY(?) AND f.date BETWEEN ? AND ? ' ...
         ' AND EXTRACT(HOUR FROM f.date) >= ? ' ...
         ' AND EXTRACT(HOUR FROM f.date) <= ? ' ...
         ' ORDER BY f.fk_id, f.date;'];
    
    filteredColumnTypes = {'int', 'real', 'real', 'date'};
    
    jSensorIds = javaArray('java.lang.Integer', length(sensorIdList));
        
    for iSensor = 1:length(sensorIdList)
        jSensorIds(iSensor) = sensorIdList(iSensor);
    end
        
    filteredSetData = ...
        {jSensorIds, startTime, endTime, startHour, endHour};
    
    filteredSetTypes = ...
        {'aint', 'date', 'date', 'int', 'int'};
    
    filteredFields = ...
        {'pemsIds', 'density', 'fraction_working', 'date'};
    
    filteredData = convertToMatlabStruct(...
        collectDataFromDB(filteredPSName, filteredPS, filteredColumnTypes, ...
                          filteredSetData, filteredSetTypes), filteredFields);

end

function [hour,minute,second] = secondsToHoursMinutesSeconds(nbSeconds)
    hour = floor(nbSeconds/3600);
    nbSeconds2 = nbSeconds-hour*3600;
    minute = floor((nbSeconds-hour*3600)/60);
    second = nbSeconds2 - minute * 60;
end

function [year, month, day, hour, minute, second] = ...
        javaTimeToHoursMinutesSeconds(jTime)
    year = jTime.getYear();
    month = jTime.getMonth();
    day = jTime.getDayOfMonth();
    hour = jTime.getHour();
    minute = jTime.getMinute();
    second = jTime.getSecond();
end

function outputStruct = shrinkStruct(inputStruct, reductionList)
    
    if isfield(inputStruct, 'pemsIds')
        outputStruct = struct();
        
        indices = [];
        
        for iId = 1:length(reductionList)
            indices = [indices; ...
                       find(inputStruct.pemsIds == ...
                            reductionList(iId))];
        end

        inputStructFields = fieldnames(inputStruct);
        
        for iField = 1:length(inputStructFields)
            fieldValue = getfield(inputStruct, inputStructFields{iField});
            outputStruct = setfield(...
                outputStruct, ...
                inputStructFields{iField}, ...
                fieldValue(indices));
        end

    else
        display('No "pemsIds" field in this struct');
    end    
end

%% Orders the struct fields according to "orderedList"
% orderedList is a list of links.
function outputStruct = orderStruct(inputStruct, orderedList)

    if isfield(inputStruct, 'linkIds')
        outputStruct = struct();
        
        indices = zeros(length(inputStruct.linkIds), 1);
        
        structCount = 1;
        for iLink = 1:length(orderedList)
            index = find(inputStruct.linkIds == ...
                         orderedList(iLink)); 
            
            if isfield(inputStruct, 'offset') & length(index) > 1
                
                sortedOffsets = sort(inputStruct.offset(index));
                
                for iIndex = 1:length(index)
                    indices(structCount) = ...
                        find(inputStruct.offset ...
                             == sortedOffsets(iIndex));
                    structCount = structCount + 1;
                end
                
            else
                
                for iIndex = 1:length(index)
                
                    indices(structCount) = index(iIndex);
                    structCount = structCount + 1;
                end
            end
        end
        
        inputStructFields = fieldnames(inputStruct);
        
        for iField = 1:length(inputStructFields)
            fieldValue = getfield(inputStruct, inputStructFields{iField});
            outputStruct = setfield(...
                outputStruct, ...
                inputStructFields{iField}, ...
                fieldValue(indices));
        end
    
    else
        display('No "linkIds" fields in the struct');
    end
    
end